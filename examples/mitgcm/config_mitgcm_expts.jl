
import PALEOboxes as PB

import PALEOocean
import PALEOmodel
using Printf

include("../atmreservoirreaction.jl") # temporary solution to make ReactionReservoirAtm available


function build_outfilename(outfileroot, segment_number)
    return outfileroot*@sprintf("_%04i", segment_number)
end

function concatenate_segments(outfileroot, segment_range)
    segment_range = collect(segment_range)

    output = PALEOmodel.OutputWriters.load_jld2!(PALEOmodel.OutputWriters.OutputMemory(),
                    build_outfilename(outfileroot, first(segment_range)))
    for i in 2:length(segment_range)
        segment_number = segment_range[i]
        segout = PALEOmodel.OutputWriters.load_jld2!(PALEOmodel.OutputWriters.OutputMemory(),
                        build_outfilename(outfileroot, segment_number))
        
        append!(output, segout)
      
    end

    return output
end

