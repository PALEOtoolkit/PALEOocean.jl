module PALEOocean

import SnoopPrecompile

import Logging
import PALEOboxes as PB

function moduledir()
    return dirname(@__DIR__)
end

include("ocean/Ocean.jl")

include("oceansurface/Oceansurface.jl")


end # module PALEOocean
