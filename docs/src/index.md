# PALEOocean.jl documentation

## Installation and running the model

### Installation

NB: requires Julia 1.6 or later.  To check the Julia version:

    julia> versioninfo()

Clone this github repository to local directory `PALEOocean`: from a linux bash prompt or a Windows terminal,

    $ git clone https://github.com/PALEOtoolkit/PALEOocean.jl.git PALEOcopse

Start julia and navigate to the `PALEOocean/examples` folder, and run `setup.jl` to configure the `PALEOocean/examples`
Julia environment to use the local (downloaded) version of the PALEOocean package:

    julia> cd("PALEOocean/examples")
    julia> include("setup.jl") # use the local version of PALEOocean packages to allow local modifications
   
### Running a minimal transport example
Start julia and navigate to the `PALEOocean` folder, then:

    julia> cd("examples/transport_examples")
    julia> import Pkg
    julia> Pkg.activate("..") # use the PALEOocean/examples environment

    julia> include("PALEO_examples_transport_advect.jl")


## Using PALEOocean Reactions from other models

The PALEO Reactions comprising the PALEOocean models are available when the registered PALEOocean package is loaded (without downloading the repository), ie

    julia> Pkg.add("PALEOocean")
    julia> import PALEOocean


