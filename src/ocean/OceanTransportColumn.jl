module OceanTransportColumn

import LinearAlgebra
import NCDatasets

# import Infiltrator
import PALEOboxes as PB
using PALEOboxes.DocStrings

import ..Ocean

"""
    ReactionOceanTransportColumn

Set up 1D ocean column grid, provide tracer transport using supplied eddy diffusivity `Kz`.

A netCDF file (name specified by `grid_file` parameter) should provide 
`zupper`, `zmid`, `zlower` (m) defining z coordinates of cell surfaces and centres for a 1D column.

Eddy diffusivity on cell upper surfaces should be provided by a variable `Kz` (m^2 s-1).

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_SETUP)
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionOceanTransportColumn{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(        

        PB.ParString("grid_file", "S2P3_depth80_m2amp04.nc",
            description="netcdf file with grid data (zmid, zupper, zlower)"),

        PB.ParDouble("column_area", 1.0, units="m^2",
            description="column area"),
    )

 
    grid_ocean          = nothing
    grid_oceansurface   = nothing
    grid_oceanfloor     = nothing

end
   


function PB.register_methods!(rj::ReactionOceanTransportColumn, model::PB.Model)

    PB.add_method_setup!(rj, setup_grid, (PB.VarList_namedtuple(Ocean.grid_vars_all),))

    return nothing
end

function PB.set_model_geometry(rj::ReactionOceanTransportColumn, model::PB.Model)

    NCDatasets.NCDataset(rj.pars.grid_file[], "r") do nc_phys
   
        ocean_cells = nc_phys.dim["cells"]
        @info "set_model_geometry: $(PB.fullname(rj)) defining 1 column with $ocean_cells cells from grid information in file $(rj.pars.grid_file[])"

        rj.grid_ocean = PB.Grids.UnstructuredColumnGrid(                       
            ncells=ocean_cells,
            columnnames=[:-],
            Icolumns=[collect(1:ocean_cells)]
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
            cellnames=Dict(rj.grid_ocean.columnnames[i]=>i for i in 1:length(rj.grid_ocean.columnnames)))
        PB.Grids.set_subdomain!(rj.grid_oceansurface, "ocean", PB.Grids.InteriorSubdomain(ocean_cells, isurf), true)

        rj.grid_oceanfloor = PB.Grids.UnstructuredVectorGrid(
            ncells=length(ifloor),
            cellnames=Dict(rj.grid_ocean.columnnames[i]=>i for i in 1:length(rj.grid_ocean.columnnames)))
        PB.Grids.set_subdomain!(rj.grid_oceanfloor, "ocean", PB.Grids.InteriorSubdomain(ocean_cells, ifloor), true)

        Ocean.set_model_domains(model, rj.grid_ocean, rj.grid_oceansurface, rj.grid_oceanfloor)
    end   
    
    return nothing
end

function setup_grid(
    m::PB.ReactionMethod,
    pars,
    (grid_vars, ),
    cellrange::PB.AbstractCellRange,
    attribute_name
)
    attribute_name == :setup || return

    rj = m.reaction

    @info "setup_grid: $(PB.fullname(rj)) reading grid from file $(pars.grid_file[])"

    NCDatasets.NCDataset(pars.grid_file[], "r") do nc_phys
   
        grid_vars.zupper        .= nc_phys["zupper"][:]
        grid_vars.zlower        .= nc_phys["zlower"][:]
        grid_vars.zmid          .= nc_phys["zmid"][:]

        # attach coordinates to grid for output visualisation etc
        if isdefined(PB, :set_coordinates!) # PALEOboxes >= 0.22
            PB.set_coordinates!(rj.domain.grid, "cells", ["zmid", "zlower", "zupper"])
        else
            empty!(rj.domain.grid.z_coords)
            push!(rj.domain.grid.z_coords, PB.FixedCoord("zmid", grid_vars.zmid, PB.get_variable(m, "zmid").attributes))
            push!(rj.domain.grid.z_coords, PB.FixedCoord("zlower", grid_vars.zlower, PB.get_variable(m, "zlower").attributes))
            push!(rj.domain.grid.z_coords, PB.FixedCoord("zupper", grid_vars.zupper, PB.get_variable(m, "zupper").attributes))
        end

        grid_vars.volume        .= pars.column_area[].*(grid_vars.zupper .- grid_vars.zlower)
        grid_vars.volume_total[]= sum(grid_vars.volume)

        grid_vars.Abox          .= pars.column_area[]
        grid_vars.Asurf[1]      = pars.column_area[]

        grid_vars.Afloor[1]     = pars.column_area[]
        grid_vars.Afloor_total[]= sum(grid_vars.Afloor)
        grid_vars.zfloor[1]     = grid_vars.zlower[end]
        
        # constant density
        grid_vars.rho_ref       .= 1027

        grid_vars.pressure      .= -grid_vars.zmid   # pressure(dbar) ~ depth (m)
    end

    return nothing
end     

           
function PB.register_dynamic_methods!(rj::ReactionOceanTransportColumn)

    (transport_conc_vars, transport_sms_vars, _) =
        Ocean.find_transport_vars(rj.domain, add_transport_input_vars=false)

    vars = [
        PB.VarDep("Kz",    "m^2 s-1", "Vertical turbulent diffusivity. NB: defined on upper cell faces"),
    ]
 
    PB.add_method_do!(
        rj, 
        do_ocean_1D_column_trspt,
        (   
            PB.VarList_namedtuple(vars),
            PB.VarList_namedtuple(PB.VarDep.(Ocean.grid_vars_ocean)),
            PB.VarList_components(transport_conc_vars),
            PB.VarList_components(transport_sms_vars),            
        );
        preparefn=prepare_do_ocean_1D_column_trspt,
    )

    return nothing
end

# add a buffer to store transport matrix, of appropriate element type (which may be an AD type)
function prepare_do_ocean_1D_column_trspt(
    m::PB.ReactionMethod, (vars, grid_vars, transport_conc_components, transport_sms_components)
)
    rj = m.reaction

    ET = eltype(vars.Kz)

    dtm_buffer = LinearAlgebra.Tridiagonal(
        zeros(ET, rj.grid_ocean.ncells-1), zeros(ET, rj.grid_ocean.ncells), zeros(ET, rj.grid_ocean.ncells-1)
    )

    return (vars, grid_vars, transport_conc_components, transport_sms_components, dtm_buffer)
end

function do_ocean_1D_column_trspt(
    m::PB.ReactionMethod,
    pars,
    (vars, grid_vars, transport_conc_components, transport_sms_components, dtm_buffer),
    cellrange::PB.AbstractCellRange,
    deltat
)

    dtm_buffer.d .= 0
    for i in cellrange.indices
        # calculate transport between cell i and i+1
        # NB: Kz is defined on upper cell faces
        i == last(cellrange.indices) && continue

        #   m^3 s-1      m^2 s-1            m^2                  / m
        kdtm    = vars.Kz[i+1]*pars.column_area[]/(grid_vars.zmid[i]-grid_vars.zmid[i+1])
        #  yr-1         m^3 s-1  /   m^3              * s yr-1
        dtm_buffer.dl[i] = kdtm/grid_vars.volume[i+1]*PB.Constants.k_secpyr
        dtm_buffer.du[i] = kdtm/grid_vars.volume[i]*PB.Constants.k_secpyr
        dtm_buffer.d[i] -= dtm_buffer.du[i] 
        dtm_buffer.d[i+1] -= dtm_buffer.dl[i] 
    end


    Ocean.do_transport(
        grid_vars, transport_conc_components, transport_sms_components,
        dtm_buffer,
        cellrange
    )

    return nothing
end

end # module
