# PALEOocean.jl

[![CI](https://github.com/PALEOtoolkit/PALEOocean.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/PALEOtoolkit/PALEOocean.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://PALEOtoolkit.github.io/PALEOocean.jl/dev)

Ocean components for the PALEO biogeochemical model. 


**NB: work-in-progress - this repo contains initial minimal examples only to test infrastructure.**

## Installation and running a minimal example

### Installation

NB: requires Julia 1.6 or later.  To check the Julia version:

    julia> versioninfo()

Clone this github repository to local directory `PALEOocean`:

from a linux bash prompt or a Windows terminal,

    $ git clone https://github.com/PALEOtoolkit/PALEOocean.jl.git PALEOocean

Start julia and navigate to the `PALEOocean/examples` folder, and run `setup.jl` to configure the `PALEOocean/examples`
Julia environment to use the local (downloaded) version of the PALEOocean package:

    julia> cd("PALEOocean/examples")
    julia> include("setup.jl") # use the local version of PALEOocean packages to allow local modifications
   
### Running a minimal ocean transport example
Start julia and navigate to the `PALEOocean` folder, then:

    julia> cd("PALEOocean/examples")
    julia> import Pkg
    julia> Pkg.activate(".") # use the PALEOocean/examples environment

    julia> cd("transport_examples")
    julia> include("PALEO_examples_transport_advect.jl")


## Using PALEOocean Reactions from other models

The PALEOocean Reactions are available to the PALEO framework when the registered PALEOocean package is loaded (without downloading the repository), ie

    julia> Pkg.add("PALEOocean")
    julia> import PALEOocean