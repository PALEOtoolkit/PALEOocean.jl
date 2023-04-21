# 3 Box Ocean Examples

These examples demonstrate the 3-box [Sarmiento1984](@cite), [Toggweiler1985](@cite) ocean model,
standalone and coupled to the COPSE land surface and sediment/crust (as used in [Clarkson2015](@cite)).

## Abiotic CO2/DIC only atmosphere-ocean 

    julia> include("PALEO_examples_oaonly_abiotic.jl")

Abiotic atmosphere-ocean with atmosphere CO2, ocean DIC, TAlk. Test case cf Sarmiento & Toggweiler (2007) book, Fig 10.4, p436-7

Commented-out options in file to set k_piston to show effect of default/fast/slow air-sea exchange rates

## Biotic atmosphere-ocean (no weathering or burial)

    julia> include("PALEO_examples_oaonly.jl")

Biotic atmosphere-ocean with atmosphere O2, CO2, ocean P, O2, SO4/H2S, CH4, DIC, TAlk (no weathering or burial).

Use in conjunction with expt='killbio' (disables production at t=0 yr) to
demonstrate effect of biological pump.

## Open atmosphere-ocean with silicate/carbonate weathering and burial

    julia> include("PALEO_examples_oaopencarb.jl")

Biotic atmosphere-ocean with atmosphere O2, CO2, ocean P, O2, SO4/H2S, CH4, DIC, TAlk

Open atm-ocean carbonate system, with carbonate/silicate weathering input, degassing input, and carbonate burial output

Closed ocean organic carbon, sulphur systems (no burial)
