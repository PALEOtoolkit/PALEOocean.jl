module OceanTransport3box

import SparseArrays
import PALEOboxes as PB
using PALEOboxes.DocStrings

import PALEOocean

# import Infiltrator # Julia debugger

"""
    ReactionOceanTransport3box

3-box [Sarmiento1984](@cite), [Toggweiler1985](@cite) ocean model.

    ---------------------------------------
    |  1 (s)                 | 2(h)        |
    |                    -->--->-          |
    |-------------------|----|  |          |   
    |                   |    |  |     /|\\  |
    |                   |    |__|______|___|    
    |                   |       |      |   |
    |  3 (d)             -<--<--      \\|/  |
    |                  circT          fhd  |
    |                                      |
    ----------------------------------------

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_SETUP)
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionOceanTransport3box{P} <: PB.AbstractReaction
    base::PB.ReactionBase
  
    pars::P = PB.ParametersTuple(
        PB.ParDouble("circT",     20*1e6, units="m^3 s^-1", 
            description="overturning circulation"),
        PB.ParDouble("circfhd",   60*1e6, units="m^3 s^-1", 
            description="high latitude <-> deep exchange rate"),
        PB.ParDoubleVec("temp",  [21.5, 2.5, 2.5], units="degrees C", 
            description="ocean temperature"),
        PB.ParBool("temp_trackglobal", false,
            description="track global temperature (apply offset of global temp -15C"),
    )

    "Ocean circulation (defined as transport matrix)
     NB: tracers are column vectors, multiplied by transport
            
    Transport matrix:  [1->  2->1  3->1
                         1->2 2->   3->2
                         1->3 2->3  3-> ]
    For internal use: Units: m^3 s^-1"
    trspt_circ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    "Transport matrix: Units yr^{-1}
     so dc/dt = trspt_dtm * c [yr^-1]
     where c is column vector of tracer concentrations"
    trspt_dtm::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    # calculate sparse transpose to provide a test case (not an optimisation we need, given this is a 3x3 matrix!)
    trspt_dtm_tr::SparseArrays.SparseMatrixCSC{Float64, Int64} = SparseArrays.spzeros(0, 0)
end
   

function PB.set_model_geometry(rj::ReactionOceanTransport3box, model::PB.Model)
    
    ocean_cells = 3 # Number of cells (= ocean Domain size)

    # Define some named cells for plotting (only)
    oceancellnames=[:s, :h, :d]

    isurf=[1,2]
    ifloor=[1, 2, 3]

    # set minimal grid (just names for plotting)
    oceangrid = PB.Grids.UnstructuredVectorGrid(ncells=ocean_cells,                                              
                                                cellnames=Dict(oceancellnames[i]=>i for i in 1:length(oceancellnames)))

    PB.Grids.set_subdomain!(oceangrid, "oceansurface", PB.Grids.BoundarySubdomain(isurf), true)
    @info "  set ocean.oceansurface Subdomain size=$(length(isurf))"

    PB.Grids.set_subdomain!(oceangrid, "oceanfloor", PB.Grids.BoundarySubdomain(ifloor), true)
    @info "  set ocean.oceanfloor Subdomain size=$(length(ifloor))"

    surfacedom_cellnames = Dict(oceancellnames[isurf[i]]=>i for i in eachindex(isurf))
    surfacegrid = PB.Grids.UnstructuredVectorGrid(ncells=length(isurf),                                             
                                                  cellnames=surfacedom_cellnames)
    PB.Grids.set_subdomain!(surfacegrid, "ocean", PB.Grids.InteriorSubdomain(ocean_cells, isurf), true)

    floordom_cellnames = Dict(oceancellnames[ifloor[i]]=>i for i in eachindex(ifloor))
    floorgrid = PB.Grids.UnstructuredVectorGrid(ncells=length(ifloor),                                              
                                                cellnames=floordom_cellnames)
    PB.Grids.set_subdomain!(floorgrid, "ocean", PB.Grids.InteriorSubdomain(ocean_cells, ifloor), true)

    PALEOocean.Ocean.set_model_domains(model, oceangrid, surfacegrid, floorgrid)    
    
    return nothing
end



function PB.register_methods!(rj::ReactionOceanTransport3box)
   
    physvars = [
        PB.VarProp("sal",                           "psu",      "Ocean salinity"),
        PB.VarProp("rho",                           "kg m^-3",  "physical ocean density"),
        # no length check as create and set oceansurface Variable from atm Domain
        PB.VarProp("oceansurface.open_area_fraction","",        "fraction of area open to atmosphere", attributes=(:check_length=>false,)),
    ]

    PB.add_method_setup!(
        rj, 
        do_setup_grid,
        (   
            PB.VarList_namedtuple(PALEOocean.Ocean.grid_vars_all), 
            PB.VarList_namedtuple(physvars),
        ),
    )

    tempvars = [
        PB.VarDepScalar("(global.TEMP)",            "K",        "global mean temperature"),
        PB.VarProp("temp",                          "K",        "Ocean temperature"),
    ]

    PB.add_method_do!(
        rj, 
        do_temperature,
        (PB.VarList_namedtuple(tempvars), ),
    )    

    return nothing
end


function do_setup_grid(m::PB.ReactionMethod, pars, (grid_vars, physvars), cellrange::PB.AbstractCellRange, attribute_name)
    attribute_name == :setup || return

    rj = m.reaction

    Atot = 3.6e14 # m^2 total ocean area
    grid_vars.Abox  .= [0.85,             0.15,       1.0       ]*Atot

    grid_vars.Asurf .= grid_vars.Abox[rj.domain.grid.subdomains["oceansurface"].indices]

    # Define Afloor for all boxes so burial can still be defined from all boxes
    # (in reality, surely small but non-zero for 'surface' boxes)
    grid_vars.Afloor            .= [0.0,             0.0, grid_vars.Abox[3]]  
    grid_vars.Afloor_total[]    = sum(grid_vars.Afloor)
    # constant density
    grid_vars.rho_ref           .= 1027
             
    masstot    = 1.3697e21 # COPSE 5_14 value of ocean mass
    grid_vars.volume_total[]    = 1.3697e21/grid_vars.rho_ref[1]        # Ocean volume m^3  - specified to keep COPSE 5_14 value of ocean mass  
    grid_vars.volume[1:2]       .= [100*grid_vars.Asurf[1], 250*grid_vars.Asurf[2]] # Volume of ocean boxes, m^3
    grid_vars.volume[3]         = grid_vars.volume_total[] - sum(grid_vars.volume[1:2])
           
    grid_vars.zupper            .= -[0      ,    0 ,         100]
    grid_vars.zlower            .= grid_vars.zupper-grid_vars.volume./grid_vars.Abox
    grid_vars.zmid              .= 0.5.*(grid_vars.zupper .+ grid_vars.zlower)
    grid_vars.zfloor            .= grid_vars.zlower

    grid_vars.pressure          .= -grid_vars.zmid   # pressure(dbar) ~ depth (m)

    # set salinity
    physvars.sal                    .= 35.0
    physvars.rho                    .= 1027
    physvars.open_area_fraction     .= 1.0
    
    # transport matrix m^3 s-1
    rj.trspt_circ  = [
        -pars.circT[]       0                                  pars.circT[]
        pars.circT[]        -(pars.circT[]+pars.circfhd[])     pars.circfhd[]
        0                   pars.circT[]+pars.circfhd[]       -(pars.circT[]+pars.circfhd[])
    ]

    # convert to yr^{-1}
    rj.trspt_dtm = PB.Constants.k_secpyr*rj.trspt_circ./repeat(grid_vars.volume, 1, length(grid_vars.volume))

    # calculate sparse transpose to provide a test case (not an optimisation, given this is a 3x3 matrix!)
    rj.trspt_dtm_tr = SparseArrays.sparse(transpose(rj.trspt_dtm))

    return nothing
end     
           

function do_temperature(m::PB.ReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)

    # Set temperature
    if pars.temp_trackglobal[]
        vars.temp      .= pars.temp .- 15.0 .+ vars.TEMP[] # temperature (K)
    else
        vars.temp      .= pars.temp .+ PB.Constants.k_CtoK # temperature (K)
    end

    return nothing
end

function PB.register_dynamic_methods!(rj::ReactionOceanTransport3box)

    (transport_conc_vars, transport_sms_vars, transport_input_vars) =
        PALEOocean.Ocean.find_transport_vars(rj.domain, add_transport_input_vars=true)

    PB.add_method_do!(
        rj, 
        do_transport,
        (   
            PB.VarList_namedtuple(PB.VarDep.(PALEOocean.Ocean.grid_vars_ocean)),
            PB.VarList_components(transport_conc_vars),
            PB.VarList_components(transport_sms_vars),
            PB.VarList_components(transport_input_vars),
        ),
        preparefn=PALEOocean.Ocean.prepare_transport
    )

    return nothing
end


function do_transport(
    m::PB.ReactionMethod,
    (grid_vars, transport_conc_components, transport_sms_components, transport_input_components, buffer), 
    cellrange::PB.AbstractCellRange, 
    deltat
)
    rj = m.reaction

    PALEOocean.Ocean.do_transport(
        grid_vars, transport_conc_components, transport_sms_components, transport_input_components, buffer,
        rj.trspt_dtm,
        cellrange
    )

    # PALEOocean.Ocean.do_transport_tr(
    #     grid_vars, transport_conc_components, transport_sms_components, transport_input_components, buffer,
    #     rj.trspt_dtm_tr, 
    #     cellrange
    # ) 
    return nothing
end


end
