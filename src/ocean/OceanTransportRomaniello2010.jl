module OceanTransportRomaniello2010

import SparseArrays
import MAT   # Matlab file access

import PALEOboxes as PB
using PALEOboxes.DocStrings
import ..Ocean

import Infiltrator # Julia debugger

"""
    ReactionOceanTransportRomanielloICBM

Modern global and Black Sea ocean transport from [Romaniello2010a](@cite) [Romaniello2010](@cite),

# Global configurations (`circname = Global_79_Box`, `circname = Global_13_Box`):

The `ocean` Domain consists of three columns, :hlat, :gyre, :upw

Ocean box indices are:

      :hlat   :gyre   :upw
    ---------------------------
    |  1    |  14   |  47     |
     |      |       |        |
      |     |       |       | 
       |    |       |      |
       | 13 |  46   | 79  |  
        -------------------

      :hlat   :gyre   :upw
    ---------------------------
    |  1    |  4    |  9      |
     |      |       |        |
      |     |       |       | 
       |    |       |      |
       | 3  |  8    | 13  |  
       -------------------

The `oceansurface` Domain has three boxes, the `oceanfloor` Domain has 1 box per ocean box.

# Black Sea configuration  (`circname = Black_Sea`)

The `ocean` Domain consists of a single column. The `oceansurface` Domain has one box, the `oceanfloor` Domain has 1 box per ocean box.

Bosphorus outflow flux (Parameter `bosph_outflow`) is represented by Variables with default linking to a `fluxBosphorusOutflow` Domain:

    bosph_outflow_<X> --> fluxBosphorusOutflow.flux_<X>

The .yaml configuration file can be used to configure a closed system by adding outflow flux back into the ocean surface:

    variable_links:
        bosph_outflow_*: ocean.oceansurface.*_sms

Bosphorus inflow (plume) flux is represented by Variables with default linking to concentrations and source - sink fluxes in a `BosphorusInflow` Domain:

    bosph_inflow_<X>_conc --> BosphorusInflow.<X>_conc
    bosph_inflow_<X>_sms --> BosphorusInflow.<X>_sms    # will be -ve for flux into Black Sea

the .yaml configuration file can be used to configure a closed system by sourcing inflow flux from ocean surface:

    variable_links:
        bosph_inflow_*_conc: ocean.oceansurface.*_conc
        bosph_inflow_*_sms: ocean.oceansurface.*_sms


# Implementation
Reads Matlab .mat files created from the published SI with Matlab commands:

    >> circ_global_79_box = Global_79_Box_ICBM_Params
    >> circ_global_13_box = Global_13_Box_ICBM_Params
    >> circ_black_sea = Black_Sea_ICBM_Params
    >> save('romaniello_global79','-struct', 'circ_global_79_box',  '-v6')
    >> save('romaniello_global13','-struct', 'circ_global_13_box',  '-v6')
    >> save('romaniello_blacksea','-struct', 'circ_black_sea',  '-v6')

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_SETUP)
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionOceanTransportRomanielloICBM{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParString("matdir", "romaniello2010_transport",
            description="folder with Romaniello (2010) transport and geometry data files"),
        PB.ParString("circname", "Global_79_Box", allowed_values=["Global_79_Box", "Global_13_Box", "Black_Sea"],
            description="transport matrix from Romaniello(2010) G3"),

        PB.ParBool("set_temp", true,
            description="true to provide and set ocean.temp variable to (very approximate) values"),

        PB.ParBool("temp_trackglobal", false,
            description="track global temperature (apply offset of global temp -15C"),

        PB.ParDouble("bosph_outflow", 6.04e11, units="m^3 yr^-1",
            description="Bosphorus outflow (Black_Sea only)"),

        PB.ParDouble("bosph_inflow", 3.05e11, units="m^3 yr^-1",
            description="Bosphorus inflow (plume, Black_Sea only)"),
    )

 
    temp_oceanC::Vector{Float64} = Float64[] # ocean temperature (degrees C)

    grid_ocean          = nothing
    grid_oceansurface   = nothing
    grid_oceanfloor     = nothing
    

    "Transport matrix: Units yr^{-1}
     so dc/dt = trspt_dtm * c [yr^-1]
     where c is column vector of tracer concentrations"
    trspt_dtm::SparseArrays.SparseMatrixCSC{Float64,Int64} = SparseArrays.spzeros(0, 0)

    "Transpose of trspt_dtm"
    trspt_dtm_tr::SparseArrays.SparseMatrixCSC{Float64,Int64} = SparseArrays.spzeros(0, 0)

    "Black_Sea only: fraction of Bosphorus inflow (plume) flux to each box"
    plume_frac::SparseArrays.SparseVector{Float64, Int64} = SparseArrays.spzeros(0)

    rom_data                    = nothing  # raw data from .mat file
end
   


function read_datafiles(rj::ReactionOceanTransportRomanielloICBM)

    matdir = rj.pars.matdir[] # directory containing .mat files

    # rom keys:
    # rom[<colname>]                        indices of boxes for this column, ordered surface to floor
    # rom["depths"]                         mid depths (m) cells
    # rom[<colname>_bnd]                    indices of horizontal surfaces for this column (length of column + 1)
    # rom["bnd_depths"][rom[<colname>_bnd]] depths (m) of horizontal surfaces (length of column + 1)
    # rom["Hyps"][rom[<colname>_bnd]]       area (m^2) of horizontal surfaces (length of column + 1)
    # rom["T"]                              transport matrix

    if rj.pars.circname[] == "Global_13_Box"
        matfilename = joinpath(matdir, "romaniello_global13.mat")
        @info "read_datafiles: $(fullname(rj)) reading transport matrix from file $(matfilename)"
        rj.rom_data = MAT.matread(matfilename)
        rj.grid_ocean = PB.Grids.UnstructuredColumnGrid(
            # domain=rj.domain,
            ncells=rj.rom_data["nboxes"],
            columnnames=[:hlat, :gyre, :upw],
            Icolumns=[Int.(vec(rj.rom_data["hlat"])), Int.(vec(rj.rom_data["gyre"])), Int.(vec(rj.rom_data["upw"]))])
    elseif rj.pars.circname[] == "Global_79_Box"
        matfilename = joinpath(matdir, "romaniello_global79.mat")
        @info "read_datafiles: $(PB.fullname(rj)) reading transport matrix from file $(matfilename)"
        rj.rom_data = MAT.matread(matfilename)
        rj.grid_ocean = PB.Grids.UnstructuredColumnGrid(
            # domain=rj.domain,
            ncells=rj.rom_data["nboxes"],
            columnnames=[:hlat, :gyre, :upw],
            Icolumns=[Int.(vec(rj.rom_data["hlat"])), Int.(vec(rj.rom_data["gyre"])), Int.(vec(rj.rom_data["upw"]))])
    elseif rj.pars.circname[] == "Black_Sea"
        matfilename = joinpath(matdir, "romaniello_blacksea.mat")
        @info "read_datafiles: $(PB.fullname(rj)) reading transport matrix from file $(matfilename)"
        rj.rom_data = MAT.matread(matfilename)
        # single column, no name
        rj.grid_ocean = PB.Grids.UnstructuredColumnGrid(
            # domain=rj.domain,
            ncells=rj.rom_data["nboxes"],
            columnnames=[:-],
            Icolumns=[collect(1:rj.rom_data["nboxes"])])
    else
        error("unknown circname '$(rj.pars.circname[])")
    end
    
    return nothing
end


function PB.set_model_geometry(rj::ReactionOceanTransportRomanielloICBM, model::PB.Model)

    read_datafiles(rj)

    ocean_cells = rj.grid_ocean.ncells

    # set subdomain mappings
    isurf=[col[1] for col in rj.grid_ocean.Icolumns]  # first index in each column is surface box
    ifloor=collect(1:ocean_cells)  # every box has an oceanfloor box under it

    PB.Grids.set_subdomain!(rj.grid_ocean, "oceansurface", PB.Grids.BoundarySubdomain(isurf), true)
    @info "  set ocean.oceansurface Subdomain size=$(length(isurf))"
    PB.Grids.set_subdomain!(rj.grid_ocean, "oceanfloor", PB.Grids.BoundarySubdomain(ifloor), true)
    @info "  set ocean.oceanfloor Subdomain size=$(length(ifloor))"

    rj.grid_oceansurface = PB.Grids.UnstructuredVectorGrid(
        # domain=PB.get_domain(model, "oceansurface"),
        ncells=length(isurf),
        cellnames=Dict(rj.grid_ocean.columnnames[i]=>i for i in 1:length(rj.grid_ocean.columnnames)))
    PB.Grids.set_subdomain!(rj.grid_oceansurface, "ocean", PB.Grids.InteriorSubdomain(ocean_cells, isurf), true)

    rj.grid_oceanfloor = PB.Grids.UnstructuredColumnGrid(
        # domain=PB.get_domain(model, "oceanfloor"),
        ncells=length(ifloor),
        columnnames=rj.grid_ocean.columnnames,
        Icolumns=rj.grid_ocean.Icolumns)
    PB.Grids.set_subdomain!(rj.grid_oceanfloor, "ocean", PB.Grids.InteriorSubdomain(ocean_cells, ifloor), true)

    Ocean.set_model_domains(model, rj.grid_ocean, rj.grid_oceansurface, rj.grid_oceanfloor)    
    
    return nothing
end

function PB.register_methods!(rj::ReactionOceanTransportRomanielloICBM)
   
    physvars = [
        PB.VarProp("sal",                 "psu", "Ocean salinity"),
        PB.VarProp("rho", "kg m^-3", "physical ocean density"),
        # no length check as create and set oceansurface Variable from atm Domain
        PB.VarProp("oceansurface.open_area_fraction","","fraction of area open to atmosphere", attributes=(:check_length=>false,)),
    ]

    PB.add_method_setup!(
        rj, 
        do_setup_grid,
        (   
            PB.VarList_namedtuple(Ocean.grid_vars_all), 
            PB.VarList_namedtuple(physvars),
        ),
    )

    if rj.pars.set_temp[]
        tempvars = [
            PB.VarDepScalar("(global.TEMP)",            "K",        "global mean temperature"),
            PB.VarProp("temp",                          "K",        "Ocean temperature"),
        ]

        PB.add_method_do!(
            rj, 
            do_temperature,
            (PB.VarList_namedtuple(tempvars), ),
        )
    end

    return nothing
end


function do_setup_grid(
    m::PB.ReactionMethod,
    (grid_vars, physvars),
    cellrange::PB.AbstractCellRange,
    attribute_name
)
    attribute_name == :setup || return
    
    rj = m.reaction
    # rom keys:
    # rom[<colname>]                        indices of boxes for this column, ordered surface to floor
    # rom["depths"]                         mid depths (m) cells
    # rom[<colname>_bnd]                    indices of horizontal surfaces for this column (length of column + 1)
    # rom["bnd_depths"][rom[<colname>_bnd]] depths (m) of horizontal surfaces (length of column + 1)
    # rom["Hyps"][rom[<colname>_bnd]]       area (m^2) of horizontal surfaces (length of column + 1)
    # rom["T"]                              transport matrix

    for icl in 1:length(rj.grid_ocean.Icolumns)
        iboxes = rj.grid_ocean.Icolumns[icl]  # indices for this column
        if rj.pars.circname[] == "Black_Sea"
            ibnd = collect(1:length(iboxes)+1)
        else
            ibnd = Int.(rj.rom_data[String(rj.grid_ocean.columnnames[icl])*"_bnd"])
        end
    
        grid_vars.volume[iboxes] .= rj.rom_data["V"][iboxes]/1000.0 # convert to m^3

        grid_vars.Abox[iboxes]  .= rj.rom_data["Hyps"][ibnd[1:end-1]]
        grid_vars.Asurf[icl]  = rj.rom_data["Hyps"][ibnd[1]]
        # assume column closed at bottom box (Romaniello hypsometry has small term here)
        grid_vars.Afloor[iboxes[1:end-1]]  .= grid_vars.Abox[iboxes[1:end-1]] - grid_vars.Abox[iboxes[2:end]]
        grid_vars.Afloor[iboxes[end]] = grid_vars.Abox[iboxes[end]]

        # NB: there is an inconsistency in Matlab (rj.rom_data) vs PALEO conventions: volume != Abox * (zupper - zlower) 
        grid_vars.zupper[iboxes]  .= -rj.rom_data["bnd_depths"][ibnd[1:end-1]]
        grid_vars.zmid[iboxes] .= -rj.rom_data["depths"][iboxes]
        # grid_vars.zlower[iboxes]  .= grid_vars.zupper[iboxes]-grid_vars.volume[iboxes]./grid_vars.Abox[iboxes]
        grid_vars.zlower[iboxes]  .= -rj.rom_data["bnd_depths"][ibnd[2:end]]

        grid_vars.zfloor[iboxes]  .= grid_vars.zlower[iboxes]
        # m^2 (k_nbox,k_nbox) Area(i, j) of horizontal contact between lower surface of box i and box j        
        # rj.grid_vars.Aoverlap
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
    
    # read transport matrix
    rj.trspt_dtm = copy(rj.rom_data["T"])

    
    if rj.pars.circname[] == "Black_Sea"
        
        # Add Bosphorus outflow back into box 1 so trspt_dtm represents a closed system
        # (the Matlab model included this flux in the T matrix - we will explicitly add it in separately if needed)
        #      yr-1                 m^3 yr-1                m^3
        rj.trspt_dtm[1, 1] += rj.rom_data["Bosph_Outflow"]/grid_vars.volume[1]

        # Calculate the distribution of Bosphorus inflow plume to each box
        # (the Matlab model supplies this as Tplume, a transport matrix)
        # We want plume_frac, the fraction of inflow flux to each box, where sum(plume_frac) == 1.0
        # m^3 yr-1 =            m^3       *       yr-1
        plume_flux = sum(grid_vars.volume.*rj.rom_data["Tplume"][:, 1]) 
        # unitless    =       m^3       *  yr-1                    / m3 yr-1
        rj.plume_frac = grid_vars.volume.*rj.rom_data["Tplume"][:, 1] ./ plume_flux
        
        # check conservation
        # for i in 1:rj.grid_ocean.ncells
        #    println("cell ", i, " sum of flux with conc in cell i = 1 mol m-3: ", sum(rj.trspt_dtm[:, i].*grid_vars.volume))
        # end
    end
    
    rj.trspt_dtm_tr = SparseArrays.sparse(transpose(rj.trspt_dtm))

    # set salinity
    physvars.sal                .= 35.0
    # constant density
    physvars.rho                .= 1027
    
    physvars.open_area_fraction .= 1.0

     # optionally set a very approximate default ocean temperature 
     if rj.pars.set_temp[]
        rj.temp_oceanC = Vector{Float64}(undef, rj.grid_ocean.ncells)
        rj.temp_oceanC .= 2.0 # all interior boxes at 2C
        if rj.pars.circname[] == "Black_Sea"
            rj.temp_oceanC[rj.domain.grid.subdomains["oceansurface"].indices] .= 7.5 # guessed 7.5C surface temp
        else
            rj.temp_oceanC[rj.domain.grid.subdomains["oceansurface"].indices] .= [2.5, 25.0, 25.0] # guessed
        end
    end

    return nothing
end     
           


function PB.register_dynamic_methods!(rj::ReactionOceanTransportRomanielloICBM)

    (transport_conc_vars, transport_sms_vars, transport_input_vars) =
        Ocean.find_transport_vars(rj.domain, add_transport_input_vars=true)

    PB.add_method_do!(
        rj, 
        do_transport,
        (   
            PB.VarList_namedtuple(PB.VarDep.(Ocean.grid_vars_ocean)),
            PB.VarList_components(transport_conc_vars),
            PB.VarList_components(transport_sms_vars),
            PB.VarList_components(transport_input_vars),
        ),
        preparefn=Ocean.prepare_transport
    )

    if rj.pars.circname[] == "Black_Sea"
        tracernames = [v.localname[1:end-5] for v in transport_conc_vars]
        
        bosph_outflow_vars = [
            PB.VarContribScalar("bosph_outflow_$n"=>"fluxBosphorusOutflow.flux_$n", "mol yr-1", "Bosphorus outflow flux")
            for n in tracernames
        ]
        PB.add_method_do!(
            rj, 
            do_bosph_outflow,
            (   
                PB.VarList_namedtuple(PB.VarDep.(Ocean.grid_vars_ocean)),
                PB.VarList_components(transport_conc_vars),
                PB.VarList_components(transport_sms_vars),
                PB.VarList_components(bosph_outflow_vars),
            ),
        )

        bosph_inflow_conc_vars = [
            PB.VarDepScalar("bosph_inflow_$(n)_conc"=>"BosphorusInflow.$(n)_conc", "mol m-3", "Bosphorus inflow concentration")
            for n in tracernames
        ]
        bosph_inflow_sms_vars  = [
            PB.VarContribScalar("bosph_inflow_$(n)_sms"=>"BosphorusInflow.$(n)_sms", "mol yr-1", "Bosphorus inflow source - sink (-ve for flux to Black Sea)")
            for n in tracernames
        ]
        PB.add_method_do!(
            rj, 
            do_bosph_inflow,
            (   
                PB.VarList_components(transport_sms_vars),
                PB.VarList_components(bosph_inflow_conc_vars),
                PB.VarList_components(bosph_inflow_sms_vars),
            ),
        )
    end
    
    return nothing
end


function do_transport(
    m::PB.ReactionMethod,
    (grid_vars, transport_conc_components, transport_sms_components, transport_input_components, buffer), 
    cellrange::PB.AbstractCellRange, 
    deltat
)
    rj = m.reaction

    Ocean.do_transport_tr(
        grid_vars, transport_conc_components, transport_sms_components, transport_input_components, buffer,
        rj.trspt_dtm_tr, 
        cellrange) 
    return nothing
end

function do_bosph_outflow(
    m::PB.ReactionMethod,
    pars,
    (grid_vars, transport_conc_components, transport_sms_components, bosph_outflow_components), 
    cellrange::PB.AbstractCellRange, 
    deltat
)

    errmsg="do_bosph_outflow: components length mismatch transport_conc_components, transport_sms_components, bosph_outflow_components (check :field_data (ScalarData, IsotopeLinear etc) match)"
    for (conc, sms, outflow) in PB.IteratorUtils.zipstrict(
                                    transport_conc_components, transport_sms_components, bosph_outflow_components; errmsg=errmsg
                                )
        for i in cellrange.indices
            if i == 1 # surface box            
                # mol yr-1 =   m^3 yr-1       * mol m-3
                flux = pars.bosph_outflow[]*conc[1]
                sms[1] -= flux
                outflow[] += flux
            end
        end
    end

    return nothing
end

function do_bosph_inflow(
    m::PB.ReactionMethod,
    pars,
    (transport_sms_components, bosph_inflow_conc_components, bosph_inflow_sms_components), 
    cellrange::PB.AbstractCellRange, 
    deltat
)
    rj = m.reaction

    errmsg = "do_bosph_inflow: components length mismatch transport_sms_components, bosph_inflow_conc_components, bosph_inflow_sms_components (check :field_data (ScalarData, IsotopeLinear etc) match)"
    for (sms, inflow_conc, inflow_sms) in PB.IteratorUtils.zipstrict(
                                            transport_sms_components, bosph_inflow_conc_components, bosph_inflow_sms_components; errmsg=errmsg
                                        )
        for i in cellrange.indices                
            # mol yr-1 =   m^3 yr-1       * mol m-3
            flux = pars.bosph_inflow[]*inflow_conc[]*rj.plume_frac[i]
            sms[i] += flux
            inflow_sms[] -= flux
        end        
    end

    return nothing
end

function do_temperature(m::PB.ReactionMethod, pars, (tempvars, ), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction

    # Set temperature
    if pars.temp_trackglobal[]
        tempvars.temp      .= rj.temp_oceanC .- 15.0 .+ tempvars.TEMP[] # temperature (K)
    else
        tempvars.temp      .= rj.temp_oceanC .+ PB.Constants.k_CtoK # temperature (K)
    end

    return nothing
end

end # module
