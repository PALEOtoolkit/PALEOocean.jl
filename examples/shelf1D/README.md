# 1D shelf sea examples and test cases

## Installation and configuration
These examples require netcdf files defining an annual cycle of physical variables: 
- water column eddy diffusivity Kz, temperature, density
- surface windspeed

The example configurations assume a zip file with an illustrative 80m deep seasonally-stratifying shelf
(50N, Celtic sea) derived from the S2P3 model has been unpacked to subfolder `S2P3_transport_20240614`

(this is collated output from the Windows 'Physics Biology Model' (s2p3.exe), http://pcwww.liv.ac.uk/~jons/model.htm
see 'Introduction to the Physical and Biological Oceanography of Shelf Seas', Simpson & Sharples (2012), CUP)

## 1D shelf examples

### O2 and passive tracer test case

    julia> include("PALEO_examples_shelf1D_O2_only.jl")

Test case demonstrating O2 air-sea exchange, and mixing of fast and slow sinking passive tracers.

### Minimal phytoplankton P, O2

    julia> include("PALEO_examples_shelf1D_P_O2.jl")

Minimal single-nutrient (parameterized as phosphorus) and light-limited phytoplankton population.

### Sulphur and carbonate system

    julia> include("PALEO_examples_shelf1D_P_O2_S_CH4_carb.jl")

Oxygenic and anoxygenic phytoplankton populations, with sulphur + methane and marine carbonate system.

## 1D shelf + sediment examples

### Minimal phytoplankton water-column sediment

    julia> include("PALEO_examples_shelf1Dsed.jl")

Coupled water-column - sediment configuration, with single phytoplankton population, sulphur + methane.
