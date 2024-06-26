# PALEOocean.jl

[![CI](https://github.com/PALEOtoolkit/PALEOocean.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/PALEOtoolkit/PALEOocean.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://PALEOtoolkit.github.io/PALEOocean.jl/dev)

The PALEOocean Julia package provides:
- a catalog of ocean circulation models, including low-dimensional box models and GCM transport-matrix representations of the global ocean at resolutions up to 1 deg.
- air-sea exchange
- vertical tracer and light transport
- biological production
- parameterisations of burial fluxes

It can be used to create standalone ocean models, or combined with other components in the [PALEOtoolkit](https://github.com/PALEOtoolkit) biogeochemical model framework to create coupled Earth and exo-Earth system models.

## Documentation

Documentation is available online at https://paleotoolkit.github.io/PALEOocean.jl/

## Installation

### Using PALEOocean Reactions from other models

The PALEOocean Reactions are available to the [PALEOtoolkit](https://github.com/PALEOtoolkit) framework when the registered PALEOocean package is installed and loaded:

    julia> Pkg.add("PALEOocean") #  install PALEOocean package in the currently active Julia environment
    julia> import PALEOocean

### Running PALEOocean examples

To install and run the PALEOocean examples, clone this github repository to local directory `PALEOocean` and run the examples from the Julia REPL.

Quickstart assuming a recent Julia installation: from a linux bash prompt or a Windows terminal,

    $ git clone https://github.com/PALEOtoolkit/PALEOocean.jl.git PALEOocean

Start julia and navigate to the `PALEOocean/examples` folder, and run `setup.jl` to configure the `PALEOocean/examples`
Julia environment to use the local (downloaded) version of the PALEOocean package:

    julia> cd("PALEOocean/examples")
    julia> include("setup.jl") # use the local version of PALEOocean packages to allow local modifications
   
Examples are in subfolders of `PALEOocean/examples/` and use the `PALEOocean/examples` Julia environment.

See [Installation and getting started](https://paleotoolkit.github.io/PALEOtutorials.jl/dev/ExampleInstallConfig/)
in the [PALEOtutorials](https://github.com/PALEOtoolkit/PALEOtutorials.jl) repository for more details including installation and configuration of Julia.
