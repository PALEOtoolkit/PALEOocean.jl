# Intermediate-complexity global ocean

Intermediate-complexity global ocean examples and test cases,
using the 3-column box model for the global ocean circulation from
[Romaniello2010](@cite).

## Installation

These examples use `ReactionOceanTransportRomanielloICBM` to read Matlab data files with the
3-column box model circulation from [Romaniello2010](@cite).

The Matlab datafiles are available as a zip file from <https://github.com/PALEOtoolkit/PALEOocean.jl/releases>,
generated from the Matlab model code available as Supplementary Information to [Romaniello2010](@cite).

The examples assume the zip file has been downloaded and unpacked to subfolder `romaniello2010_transport`

## O2 only air-sea exchange and transport test

    include("PALEO_examples_romglb_O2_only.jl")

## P, O2 with parameterized export production

    include("PALEO_examples_romglb_P_O2.jl")

Phosphorus and oxygen, with a parameterization of
export production based on light and nutrient availability.

No burial fluxes (ie oceanfloor phosphorus flux is recycled into water column),
so the ocean is effectively a closed system for phosphorus.

## P, O2 with organic carbon and phosphorus burial

    include("PALEO_examples_romglb_P_O2_open.jl")

Phosphorus and oxygen, with a parameterization of
export production based on light and nutrient availability.

Burial efficiency parameterization for burial fluxes of organic carbon and phosphorus,
with ocean phosphorus restored to modern level.

## P, O2, S, DIC with organic carbon, phosphorus, and carbonate burial

    include("PALEO_examples_romglb_P_O2_S_Carb_open.jl")

Phosphorus, oxygen, sulphur and DIC, with a parameterization of
export production based on light and nutrient availability.

Burial efficiency parameterization for burial fluxes of organic carbon and phosphorus,
with ocean phosphorus restored to modern level.

Burial efficiency parameterisation of carbonate burial, with shelf carbonate burial
based on saturation state and shelf area, and deep carbonate burial based on oceanfloor
flux and saturation state. Ocean TAlk and atmosphere pCO2 are restored to constant values
(so restoring fluxes will balance CaCO3 burial).

No pyrite or gypsum burial (so ocean is a closed system for sulphur).