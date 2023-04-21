module OceanTransport6box

import SparseArrays
import PALEOboxes as PB
using PALEOboxes.DocStrings
import PALEOocean

# import Infiltrator  # Julia debugger

"""
    ReactionOceanTransport6box

4+2-box ocean model from [Daines2016](@cite)
This is based on:
- the four-box model of [Hotinski2000](@cite)
- the five-box model of [Watson1995](@cite)
- the upwelling region representation of [Canfield2006](@cite)

There are two main ingredients here:
- An open-ocean 'intermediate/thermocline' box(i) from [Hotinski2000](@cite), coupled to the low-latitude surface box via an Ekman pumping term and to high latitude box (h) via an overturning circulation term. This is a more realistic representation of the open ocean than the 3-box [Sarmiento1984](@cite), [Toggweiler1985](@cite) model.
- A two-box shelf/slope (boxes r, rc), which can be configured as 
    - an upwelling region (k_slopetype='OMZ'), cf Canfield (2006) with upwelling from the intermediate/thermocline box
    - a low-latitude shelf (k_slopetype='shelf') with exchange terms to low-latitude surface (s) and thermocline (i) boxes

    See Oxygen oases\\Box Model 20150922\\HotinskiCalc5Box.m, HotinskiConstants5Box.m, HotinskiCalcCirc6.m

           --------------------------------------------------
          | 5(r)      |  1 (s)                 | 2(h)        |
          |           |                        |             |
          |-----------|------------------------|             |   
          | 6(rc)     |                        |             |
          |           |  3 (i)                 |_____________|    
           -----------|                        |             |
             |        |------------------------              |
             |           4(d)                                |
             |                                               |
              -----------------------------------------------    

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_SETUP)
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionOceanTransport6box{P} <: PB.AbstractReaction
    base::PB.ReactionBase
    
    pars::P = PB.ParametersTuple(
        PB.ParString("slopetype",     "OMZ",  allowed_values=["OMZ", "shelf"],
            description="type of shelf circulation (boxes r, rc)"),

        PB.ParDouble("circT",     19.3e6, units="m^3 s^-1", 
            description="overturning circulation (exchange high lat (h) to deep (d) and thermocline (i)"),
        PB.ParDouble("circfhd",   48.7e6, units="m^3 s^-1",
            description="high latitude <-> deep exchange rate"),
        PB.ParDouble("circR",     20e6, units="m^3 s^-1",
            description="upwelling thermocline (i) to slope (rc) to upwell (r) to low lat surface (s)"),
        PB.ParDouble("totalEkman",60e6, units="m^3 s^-1",
            description="total wind-driven Ekman pumping"),
        PB.ParDouble("circS",     1e6, units="m^3 s^-1",
            description="upwell(r) <-> low lat surf box (s) exchange"),

        PB.ParDoubleVec("temp",  [21.5, 2.5, 2.5, 2.5, 21.5, 2.5], units="degrees C",
            description="ocean temperature"),

        PB.ParBool("temp_trackglobal", false,
            description="track global temperature (apply offset of global temp -15C"),
    )

 
    "Ocean circulation (defined as transport matrix)
     NB: tracers are column vectors, multiplied by transport
            
    Transport matrix:  [1->  2->1  3->1, 4->1, 5->1, 6-> 1
                         1->2 2->   3->2 ...
                         1->3 2->3  3-> ...
                         ...
                         1->6 ...                   6-> ]
    For internal use: Units: m^3 s^-1"
    trspt_circ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    "Transport matrix: Units yr^{-1}
     so dc/dt = trspt_dtm * c [yr^-1]
     where c is column vector of tracer concentrations"
    trspt_dtm::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    # calculate sparse transpose to provide a test case (not an optimisation, given this is a small matrix!)
    trspt_dtm_tr::SparseArrays.SparseMatrixCSC{Float64, Int64} = SparseArrays.spzeros(0, 0)
end
   

function PB.set_model_geometry(rj::ReactionOceanTransport6box, model::PB.Model)

    ocean_cells = 6 # Number of cells (= ocean Domain size)

    # Define some named cells for plotting (only)
    oceancellnames=[:s, :h, :i, :d, :r, :rc]

    isurf=[1,2, 5]
    ifloor=[1, 2, 3, 4, 5, 6]

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

function PB.register_methods!(rj::ReactionOceanTransport6box)
   
    physvars = [
        PB.VarProp("sal",                           "psu",      "Ocean salinity"),
        PB.VarProp("rho",                           "kg m^-3",  "physical ocean density"),
        # no length check as create and set oceansurface Variable from atm Domain
        PB.VarProp("oceansurface.open_area_fraction","",        "fraction of area open to atmosphere", attributes=(:check_length=>false,)),
    ]

    PB.add_method_setup!(
        rj, 
        do_setup_grid,
        (   PB.VarList_namedtuple(PALEOocean.Ocean.grid_vars_all), 
            PB.VarList_namedtuple(physvars),
        ),
    )

    tempvars = [
        PB.VarDepScalar("(global.TEMP)",            "K",        "global mean temperature"),
        PB.VarProp("temp",                          "Kelvin",   "Ocean temperature"),
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
    #                               :s              :h          :i          :d          :r          :rc
    grid_vars.Abox  .=              [0.80,           0.15,       0.80,       1.0,        0.05,       0.05       ]*Atot

    grid_vars.Asurf .= grid_vars.Abox[rj.domain.grid.subdomains["oceansurface"].indices]

    # Define Afloor for all boxes so burial can still be defined from all boxes
    # (in reality, surely small but non-zero for 'surface' boxes)
    grid_vars.Afloor            .= [0.0,             0.0,       0.0,    grid_vars.Abox[4], 0.0,     0.0]
    grid_vars.Afloor_total[]    = sum(grid_vars.Afloor)                           
    
    # constant density
    grid_vars.rho_ref           .= 1027
    
         
    masstot    = 1.3697e21 # COPSE 5_14 value of ocean mass
    grid_vars.volume_total[]    = 1.3697e21/grid_vars.rho_ref[1]        # Ocean volume m^3  - specified to keep COPSE 5_14 value of ocean mass 
    
    if pars.slopetype[] == "shelf"
        # TODO this was not used in Daines & Lenton (2016)
        error("slopetype 'shelf' not implemented")
    elseif pars.slopetype[] == "OMZ"
         # OMZ upwelling region, cf [Canfield2006](@cite)

        #                               :s              :h          :i          :d          :r          :rc
        grid_vars.zupper        .=      -[0,            0,          100,        900,        0,          100]
        grid_vars.zlower        .=      -[100,          250,       1000,        900,        100,       1000] # box 4 :d filled in later
        
        grid_vars.volume[4] = 0.0
        grid_vars.volume    .= (grid_vars.zupper .- grid_vars.zlower).*grid_vars.Abox       # box 4 temporarily zero
        grid_vars.volume[4] = grid_vars.volume_total[] - sum(grid_vars.volume)              # update volume of box 4 to get correct total
        grid_vars.zlower[4] = grid_vars.zupper[4] - grid_vars.volume[4]/grid_vars.Abox[4]   # define depth of box 4 to get correct volume
       
        grid_vars.zmid              .= 0.5.*(grid_vars.zupper .+ grid_vars.zlower)

        grid_vars.zfloor           .= grid_vars.zlower

        circMR   = 0.0
        circTR   = pars.circR[]
    else
        error("unknown slopetype ", pars.slopetype[])
    end

    grid_vars.pressure          .= -grid_vars.zmid   # pressure(dbar) ~ depth (m)

    # set salinity
    physvars.sal                    .= 35.0
    physvars.rho                    .= 1027
    physvars.open_area_fraction     .= 1.0
    
    # calculate transport matrix (m^3 s-1)
   
    #  low lat Ekman =  total Ekman pumping - margin upwelling
    circU    = pars.totalEkman[] - pars.circR[]
    # define short names for convenience
    T, U, S, R, Fhd, TR, MR     = pars.circT[], circU, rj.pars.circS[], pars.circR[], pars.circfhd[], circTR, circMR

    rj.trspt_circ   =  [
        -(U+S+R)    0               U           0           S+R         0
        0           -(2*T+Fhd)      T           Fhd+T       0           0
        U+R         T               -(U+T+R+TR) 0           0           TR
        0           Fhd+T           0           -(Fhd+T)    0           0
        S           0               0           0           -(S+R+MR)   R+MR
        0           0               R+TR        0           MR          -(R+MR+TR)
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

function PB.register_dynamic_methods!(rj::ReactionOceanTransport6box)

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

    PALEOocean.Ocean.do_transport_tr(
        grid_vars, transport_conc_components, transport_sms_components, transport_input_components, buffer,
        rj.trspt_dtm_tr, 
        cellrange)
        
    return nothing
end


end
