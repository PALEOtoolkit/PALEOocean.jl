# Black Sea intermediate-complexity box model

Black sea examples and test cases,
using the 1-column model for the Black sea circulation from
[Romaniello2010a](@cite).

## Installation

These examples use `ReactionOceanTransportRomanielloICBM` to read Matlab data files with the
Black sea circulation from [Romaniello2010a](@cite).

The Matlab datafiles are available as a zip file from <https://github.com/PALEOtoolkit/PALEOocean.jl/releases>,
generated from the Matlab model code available as Supplementary Information to [Romaniello2010](@cite).

The examples assume the zip file has been downloaded and unpacked to subfolder `../romglb/romaniello2010_transport`
(NB: a subfolder of `romglb`, not this folder `blacksea`, as the zip file contains both global and Black sea data files).
The script `download_data_files.jl` provides a function to do this:

    include("../romglob/download_romaniello2010_files.jl")

    download_romaniello2010_files()  # download and unzip



## O2 only air-sea exchange and transport test

    include("PALEO_examples_blacksea_O2_only.jl")

## P, O2 with parameterized export production

    include("PALEO_examples_blacksea_P_O2.jl")

Phosphorus and oxygen, with a parameterization of export production based on light and nutrient availability.

Configured to create a closed circulation by returning Bosphorus outflow back to surface box, 
and sourcing Bosphorus inflow from surface box.

No burial fluxes (ie oceanfloor phosphorus flux is recycled into water column),
so the Black sea is effectively a closed system for phosphorus.


## P, O2, S with parameterized export production

    include("PALEO_examples_blacksea_P_O2_SO4.jl")

Phosphorus, oxygen, sulphur, with a parameterization of
export production based on light and nutrient availability.

