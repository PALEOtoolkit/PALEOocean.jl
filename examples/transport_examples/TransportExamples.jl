module TransportExamples

import SparseArrays
import PALEOboxes as PB
using PALEOboxes.DocStrings

import PALEOocean

# import Infiltrator # Julia debugger


"""
    ReactionTransportAdvectExample

Advection around two columns

# Parameters
$(PARS)
"""
Base.@kwdef mutable struct ReactionTransportAdvectExample{P} <: PB.AbstractReaction
    base::PB.ReactionBase
  
    pars::P = PB.ParametersTuple(
        PB.ParDouble("T",     20*1e6, units="m^3 s^-1", 
            description="advective flux"),
    )

    "Transport matrix: Units yr^{-1}
    so dc/dt = trspt_dtm * c [yr^-1]
    where c is column vector of tracer concentrations
    NB: we store the transpose for computational efficiency"
   trspt_dtm_tr::SparseArrays.SparseMatrixCSC{Float64, Int64} = SparseArrays.spzeros(0, 0)

end
   

function PB.register_methods!(rj::ReactionTransportAdvectExample)
   
    PB.add_method_setup!(
        rj, 
        do_setup_advect,
        (   
            PB.VarList_namedtuple(PB.VarDep.(PALEOocean.grid_vars_ocean)),
        ),
    )

    return nothing
end


function do_setup_advect(m::PB.ReactionMethod, pars, (grid_vars, ), cellrange::PB.AbstractCellRange, attribute_name)

    rj = m.reaction

    attribute_name == :setup || return  

    isa(rj.domain.grid, PB.Grids.UnstructuredColumnGrid) ||
        error("do_setup_advect $(PB.fullname(rj)) domain.grid is not a Grids.UnstructuredColumnGrid")
    length(rj.domain.grid.Icolumns) >= 2 ||
        error("do_setup_advect $(PB.fullname(rj)) domain.grid only has 1 column")

    # calculate transport matrix
    A = zeros(rj.domain.grid.ncells, rj.domain.grid.ncells)

    # advection in a loop down column 1, up column 2
    loopindices = vcat(rj.domain.grid.Icolumns[1], reverse(rj.domain.grid.Icolumns[2]), [first(rj.domain.grid.Icolumns[1])])
    @info "do_setup_advect: $(PB.fullname(rj)) adding advective flux $(pars.T[]) (m^3 s-1) around loop $loopindices"
    PALEOocean.add_loop!(A, grid_vars.volume, pars.T[], loopindices)

    # store the transpose as a sparse matrix for computational efficiency
    #  yr-1           s yr-1                    s-1
    rj.trspt_dtm_tr = PB.Constants.k_secpyr .* SparseArrays.sparse(transpose(A))

    return nothing
end                

function PB.register_dynamic_methods!(rj::ReactionTransportAdvectExample)

    (transport_conc_vars, transport_sms_vars, transport_input_vars) =
        PALEOocean.find_transport_vars(rj.domain, add_transport_input_vars=true)

    PB.add_method_do!(
        rj, 
        do_advect,
        (   
            PB.VarList_namedtuple(PB.VarDep.(PALEOocean.grid_vars_ocean)),
            PB.VarList_components(transport_conc_vars),
            PB.VarList_components(transport_sms_vars),
            PB.VarList_components(transport_input_vars),
        ),
        preparefn=PALEOocean.prepare_transport
    )

    return nothing
end


function do_advect(
    m::PB.ReactionMethod,
    (grid_vars, transport_conc_components, transport_sms_components, transport_input_components, buffer), 
    cellrange::PB.AbstractCellRange, 
    deltat
)
    rj = m.reaction

    PALEOocean.do_transport_tr(
        grid_vars, transport_conc_components, transport_sms_components, transport_input_components, buffer,
        rj.trspt_dtm_tr,
        cellrange
    )

    return nothing
end



"""
    ReactionTransportDiffuseExample

Eddy diffusion N columns

# Parameters
$(PARS)
"""
Base.@kwdef mutable struct ReactionTransportDiffuseExample{P} <: PB.AbstractReaction
    base::PB.ReactionBase
  
    pars::P = PB.ParametersTuple(
        PB.ParDouble("Kz",     1e-5, units="m^2 s^-1", 
            description="eddy diffusivity"),
    )

    "Transport matrix: Units yr^{-1}
     so dc/dt = trspt_dtm * c [yr^-1]
     where c is column vector of tracer concentrations
     NB: we store the transpose for computational efficiency"
    trspt_dtm_tr::SparseArrays.SparseMatrixCSC{Float64, Int64} = SparseArrays.spzeros(0, 0)
end
   

function PB.register_methods!(rj::ReactionTransportDiffuseExample)
   
    PB.add_method_setup!(
        rj, 
        do_setup_diffuse,
        (   
            PB.VarList_namedtuple(PB.VarDep.(PALEOocean.grid_vars_ocean)),
        ),
    )

    return nothing
end


function do_setup_diffuse(m::PB.ReactionMethod, pars, (grid_vars, ), cellrange::PB.AbstractCellRange, attribute_name)

    rj = m.reaction

    attribute_name == :setup || return  

    isa(rj.domain.grid, PB.Grids.UnstructuredColumnGrid) ||
        error("do_setup_advect $(PB.fullname(rj)) domain.grid is not a Grids.UnstructuredColumnGrid")
    
    # calculate transport matrix
    A = zeros(rj.domain.grid.ncells, rj.domain.grid.ncells)

    # constant diffusivity
    @info "do_setup_diffuse: $(PB.fullname(rj)) adding constant diffusive flux for Kz=$(pars.Kz[]) (m^2 s-1)"
    # eddy diffusive transport calculated as the sum of 'loops' around adjacent boxes
    for (icol, icells) in enumerate(rj.domain.grid.Icolumns)  # indices for this column, order is surface -> floor
        for (icell_top, icell_bot) in PB.IteratorUtils.zipstrict(icells[begin:end-1], icells[begin+1:end])
            # Calculate volume flux F (m^3 s-1) corresponding to eddy diffusivity Kz (m^2 s-1)
            # This follows from downgradient diffusion flux density f_c of tracer conc c,
            #              f_c          = Kz ∇c
            #              mol m^-2 s-1 = m^2 s-1 * mol m-3  / m
            #              f_c          =  Kz * (c_2 - c_1) / (zmid_2 - zmid_1) 
            # then defining F_c = A * f_c  across area A,
            #              mol s-1      = m^2 * mol m^-2 s-1
            #              F_c          = A * f_c 
            #                           =         m^3 s-1              * mol m-3
            #                           = [A * Kz / (zmid_2 - zmid_1)] * (c_2 - c_1)
            # where the [quantity in brackets] is what is needed for a transport matrix.

            # overlap area (will just be Abox in this example)
            Aintf = min(grid_vars.Abox[icell_top], grid_vars.Abox[icell_bot])
            # m^3 s-1   = m^2    / m                                                   * m^2 s-1
            F           = Aintf/(grid_vars.zmid[icell_top] - grid_vars.zmid[icell_bot])*pars.Kz[]
            loopindices = (icell_top, icell_bot, icell_top)
            PALEOocean.add_loop!(A, grid_vars.volume, F, loopindices)
        end
    end

    # store the transpose as a sparse matrix for computational efficiency
    #  yr-1         =    s yr-1                 s-1    
    rj.trspt_dtm_tr = PB.Constants.k_secpyr .* SparseArrays.sparse(transpose(A))

    return nothing
end                

function PB.register_dynamic_methods!(rj::ReactionTransportDiffuseExample)

    (transport_conc_vars, transport_sms_vars, transport_input_vars) =
        PALEOocean.find_transport_vars(rj.domain, add_transport_input_vars=true)

    PB.add_method_do!(
        rj, 
        do_diffuse,
        (   
            PB.VarList_namedtuple(PB.VarDep.(PALEOocean.grid_vars_ocean)),
            PB.VarList_components(transport_conc_vars),
            PB.VarList_components(transport_sms_vars),
            PB.VarList_components(transport_input_vars),
        ),
        preparefn=PALEOocean.prepare_transport
    )

    return nothing
end


function do_diffuse(
    m::PB.ReactionMethod,
    (grid_vars, transport_conc_components, transport_sms_components, transport_input_components, buffer), 
    cellrange::PB.AbstractCellRange, 
    deltat
)
    rj = m.reaction

    PALEOocean.do_transport_tr(
        grid_vars, transport_conc_components, transport_sms_components, transport_input_components, buffer,
        rj.trspt_dtm_tr,
        cellrange
    )

    return nothing
end

#####################################################
# Set up a grid with 2 columns
#####################################################

"""
    ReactionOceanColumnGrid

N column ocean grid, each column has `area` and `depth`, with `nlayers`

N = length(`nlayers`) ie specify layers for each column.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_SETUP)
"""
Base.@kwdef mutable struct ReactionOceanColumnGrid{P} <: PB.AbstractReaction
    base::PB.ReactionBase
  
    pars::P = PB.ParametersTuple(
        PB.ParDoubleVec("area",   PB.Constants.k_SurfAreaEarth.*[0.5, 0.5], units="m^2", 
            description="column xsection  area"),
        PB.ParDouble("depth",   2000.0, units="m", 
            description="column depth"),
        PB.ParIntVec("nlayers", [10, 10], 
            description="layers per column"),  
    )

    grid_ocean = nothing
    grid_oceansurface = nothing
    grid_oceanfloor = nothing

end
   

function PB.set_model_geometry(rj::ReactionOceanColumnGrid, model::PB.Model)
    
    @info "set_model_geometry: $(PB.fullname(rj))  columns=$(length(rj.pars.nlayers)), areas=$(rj.pars.area.v), depth=$(rj.pars.depth[]), nlayers=$(rj.pars.nlayers.v)"

    length(rj.pars.area) == length(rj.pars.nlayers) || 
        error("set_model_geometry: $(PB.fullname(rj)) area and nlayers Vector Parameters are not the same length")

    ocean_cells = sum(rj.pars.nlayers) # Number of cells (= ocean Domain size)

    Icolumns, columnnames = [], []
    coltopcell, colname = 1, 'a'
    for collayers in rj.pars.nlayers
        push!(Icolumns, collect(coltopcell:coltopcell+collayers-1))
        push!(columnnames, Symbol(colname))
        coltopcell += collayers
        colname += 1
    end
    rj.grid_ocean = PB.Grids.UnstructuredColumnGrid(                       
        ncells=ocean_cells,
        columnnames=columnnames,
        Icolumns=Icolumns
    )

    # set subdomain mappings
    isurf=[col[1] for col in rj.grid_ocean.Icolumns]  # first index in each column is surface box
    ifloor=[col[end] for col in rj.grid_ocean.Icolumns]   # last index in each column is floor box

    PB.Grids.set_subdomain!(rj.grid_ocean, "oceansurface", PB.Grids.BoundarySubdomain(isurf), true)
    @info "  set ocean.oceansurface Subdomain size=$(length(isurf))"
    PB.Grids.set_subdomain!(rj.grid_ocean, "oceanfloor", PB.Grids.BoundarySubdomain(ifloor), true)
    @info "  set ocean.oceanfloor Subdomain size=$(length(ifloor))"

    rj.grid_oceansurface = PB.Grids.UnstructuredVectorGrid(
        ncells=length(isurf),
        cellnames=Dict(name=>i for (i, name) in enumerate(rj.grid_ocean.columnnames)))
    PB.Grids.set_subdomain!(rj.grid_oceansurface, "ocean", PB.Grids.InteriorSubdomain(ocean_cells, isurf), true)

    rj.grid_oceanfloor = PB.Grids.UnstructuredVectorGrid(
        ncells=length(ifloor),
        cellnames=Dict(name=>i for (i, name) in enumerate(rj.grid_ocean.columnnames)))
    PB.Grids.set_subdomain!(rj.grid_oceanfloor, "ocean", PB.Grids.InteriorSubdomain(ocean_cells, ifloor), true)

    PALEOocean.set_model_domains(model, rj.grid_ocean, rj.grid_oceansurface, rj.grid_oceanfloor)

    return nothing
end



function PB.register_methods!(rj::ReactionOceanColumnGrid)
   
    PB.add_method_setup!(
        rj, 
        do_setup_grid,
        (   
            PB.VarList_namedtuple(PALEOocean.grid_vars_all), 
        ),
    )

    return nothing
end


function do_setup_grid(m::PB.ReactionMethod, (grid_vars, ), cellrange::PB.AbstractCellRange, attribute_name)

    rj = m.reaction

    attribute_name == :setup || return nothing

    @info "do_setup_grid: $(PB.fullname(rj)) creating $(length(rj.grid_ocean.Icolumns)) ocean columns, depth=$(rj.pars.depth[]) (m), area=$(rj.pars.area.v) (m^2)"

    for (icol, icells) in enumerate(rj.grid_ocean.Icolumns)       
        dz = rj.pars.depth[]/length(icells)
        dv = rj.pars.area[icol]*dz
        # iterate over layers. order is surface -> floor
        for (ilayer, icell) in enumerate(icells)
            # ilayer is layer within this column, icell is index in all-ocean variables
            grid_vars.zupper[icell]  = -(ilayer-1)*dz
            grid_vars.zlower[icell]  = -ilayer*dz
            grid_vars.zmid[icell] = 0.5*(grid_vars.zupper[icell] + grid_vars.zlower[icell]) 
        end        
        grid_vars.volume[icells] .= dv
        grid_vars.Abox[icells]  .= rj.pars.area[icol]

        # oceansurface setup  
        grid_vars.Asurf[icol]  = rj.pars.area[icol]       
        # oceanfloor setup
        grid_vars.Afloor[icol] = rj.pars.area[icol]
        grid_vars.zfloor[icol] = grid_vars.zlower[last(icells)]
    end   
    grid_vars.volume_total[] = sum(grid_vars.volume)
    grid_vars.Afloor_total[] = sum(grid_vars.Afloor)
    
    # attach coordinates to grid for output visualisation etc
    empty!(rj.domain.grid.z_coords)
    push!(rj.domain.grid.z_coords, PB.FixedCoord("zmid", grid_vars.zmid, PB.get_variable(m, "zmid").attributes))
    push!(rj.domain.grid.z_coords, PB.FixedCoord("zlower", grid_vars.zlower, PB.get_variable(m, "zlower").attributes))
    push!(rj.domain.grid.z_coords, PB.FixedCoord("zupper", grid_vars.zupper, PB.get_variable(m, "zupper").attributes))

    # constant density
    grid_vars.rho_ref       .= 1027
    grid_vars.pressure  .= -grid_vars.zmid   # pressure(dbar) ~ depth (m)

    return nothing
end     
           
end # module