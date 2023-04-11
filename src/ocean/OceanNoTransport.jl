module OceanNoTransport

import PALEOboxes as PB
using PALEOboxes.DocStrings

import PALEOocean

# import Infiltrator # Julia debugger

"""
    ReactionOceanNoTransport

N ocean boxes (default N=1 ie a single 0-D box), each with a surface and floor, no transport.

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_SETUP)
"""
Base.@kwdef mutable struct ReactionOceanNoTransport{P} <: PB.AbstractReaction
    base::PB.ReactionBase
  
    pars::P = PB.ParametersTuple(
        PB.ParDoubleVec("area",  [1.0], units="m^2", 
            description="surface / floor area (per box)"),
        PB.ParDoubleVec("depth",  [1.0], units="m", 
            description="depth (per box)"),
    )

    
end
   

function PB.set_model_geometry(rj::ReactionOceanNoTransport, model::PB.Model)
    
    ocean_cells = length(rj.pars.area) # Number of cells (= ocean Domain size)

    length(rj.pars.depth) == ocean_cells ||
        error("ReactionOceanNoTransport $(PB.fullname(rj)) depth and area vectors must be the same length")
    
    isurf=collect(1:ocean_cells)
    ifloor=isurf

    # set minimal grid (just names for plotting)
    oceangrid = PB.Grids.UnstructuredVectorGrid(ncells=ocean_cells)

    PB.Grids.set_subdomain!(oceangrid, "oceansurface", PB.Grids.BoundarySubdomain(isurf), true)
    @info "  set ocean.oceansurface Subdomain size=$(length(isurf))"

    PB.Grids.set_subdomain!(oceangrid, "oceanfloor", PB.Grids.BoundarySubdomain(ifloor), true)
    @info "  set ocean.oceanfloor Subdomain size=$(length(ifloor))"

    surfacegrid = PB.Grids.UnstructuredVectorGrid(ncells=length(isurf))
    PB.Grids.set_subdomain!(surfacegrid, "ocean", PB.Grids.InteriorSubdomain(ocean_cells, isurf), true)

    floorgrid = PB.Grids.UnstructuredVectorGrid(ncells=length(ifloor))
    PB.Grids.set_subdomain!(floorgrid, "ocean", PB.Grids.InteriorSubdomain(ocean_cells, ifloor), true)

    PALEOocean.Ocean.set_model_domains(model, oceangrid, surfacegrid, floorgrid)    
    
    return nothing
end



function PB.register_methods!(rj::ReactionOceanNoTransport)

    PB.add_method_setup!(
        rj, 
        do_setup_grid,
        (   
            PB.VarList_namedtuple(PALEOocean.Ocean.grid_vars_all), 
        ),
    )

    return nothing
end


function do_setup_grid(m::PB.ReactionMethod, pars, (grid_vars, ), cellrange::PB.AbstractCellRange, attribute_name)
    attribute_name == :setup || return

    grid_vars.Abox              .= pars.area
    grid_vars.Asurf             .= pars.area
    grid_vars.Afloor            .= pars.area
    grid_vars.Afloor_total[]    = sum(grid_vars.Afloor)
             
    grid_vars.volume            .=  pars.area .* pars.depth
    grid_vars.volume_total[]    = sum(grid_vars.volume)
           
    grid_vars.zupper            .= 0.0
    grid_vars.zlower            .= -1.0 .*pars.depth
    grid_vars.zmid              .= 0.5.*(grid_vars.zupper .+ grid_vars.zlower)
    grid_vars.zfloor            .= -1.0 .*pars.depth

    grid_vars.pressure          .= -grid_vars.zmid   # pressure(dbar) ~ depth (m)

    # constant density
    grid_vars.rho_ref           .= 1027 # kg m-3

    return nothing
end     


end
