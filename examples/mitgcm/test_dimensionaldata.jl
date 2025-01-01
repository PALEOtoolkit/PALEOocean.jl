# Test code for DimensionalData and Makie

import DimensionalData as DD
import CairoMakie
import PALEOboxes as PB
import PALEOmodel

# generate some test output
# include("MITgcm_2deg8_PO4MMbase.jl")
# PALEOmodel.OutputWriters.save_netcdf(paleorun.output, "test.nc")

output = PALEOmodel.OutputWriters.load_netcdf!(PALEOmodel.OutputWriters.OutputMemory(), "test.nc")

# full 3D cube lon x lat x zt at last model timestep
a_P_conc = PALEOmodel.get_array(
    output, "ocean.P_conc", (tmodel=1e12, expand_cartesian=true, squeeze_all_single_dims=false);
    omit_recorddim_if_constant=false,
)

# DimensionalData
# example of Explicit with bounds array 
# https://github.com/rafaqz/DimensionalData.jl/blob/bea9b013ff10a2ec7701f6cb223a1235dda6f842/test/lookup.jl#L215-L227

# PALEO FieldArray.dims_coords[1][2] gives Vector of coordinates attached to dimension 1 (lon)
# dims_coords[1][2][1] is cell centre, dims_coords[1][2][2] is cell bounds array, size(2, n)

dim_lon = DD.X(
    DD.Lookups.Sampled(
        a_P_conc.dims_coords[1][2][1].values,
        DD.Lookups.ForwardOrdered(),
        DD.Lookups.Explicit(a_P_conc.dims_coords[1][2][2].values),
        DD.Lookups.Intervals(DD.Lookups.Start()),
        DD.Lookups.Metadata(a_P_conc.dims_coords[1][2][1].attributes),
    )
)

dim_lat = DD.Y(
    DD.Lookups.Sampled(
        a_P_conc.dims_coords[2][2][1].values,
        DD.Lookups.ForwardOrdered(),
        DD.Lookups.Explicit(a_P_conc.dims_coords[2][2][2].values),
        DD.Lookups.Intervals(DD.Lookups.Start()),
        DD.Lookups.Metadata(a_P_conc.dims_coords[2][2][1].attributes),
    )
)

dim_zt = DD.Z(
    DD.Lookups.Sampled(
        a_P_conc.dims_coords[3][2][1].values,
        DD.Lookups.ForwardOrdered(),
        DD.Lookups.Explicit(a_P_conc.dims_coords[3][2][2].values),
        DD.Lookups.Intervals(DD.Lookups.Start()),
        DD.Lookups.Metadata(a_P_conc.dims_coords[3][2][1].attributes),
    )
)

da_P_conc = DD.DimArray(a_P_conc.values, (dim_lon, dim_lat, dim_zt))

# surface plot, integer indexing for first level
fap = CairoMakie.Makie.plot(da_P_conc[Z=1])
display(fap)

# section plot, integer indexing for lon 200 deg east
fap = CairoMakie.Makie.plot(da_P_conc[X=DD.Contains(200.0)])
display(fap)

# expand z scale, to show that bounds are being correctly used
f = CairoMakie.Figure();
ax = CairoMakie.Axis(f[1, 1])
CairoMakie.ylims!(ax, -200.0, 10.0)
CairoMakie.Makie.plot!(ax, da_P_conc[X=DD.Contains(200.0)]);
display(f)
