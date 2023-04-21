
import PALEOboxes as PB

import PALEOocean
import PALEOmodel
using Printf

include("Insolation.jl")
include("../atmreservoirreaction.jl") # temporary solution to make ReactionReservoirAtm available


function config_mitgcm_expts(baseconfig, expt, extrapars=Dict())


    if baseconfig == "abiotic_O2"

        model = PB.create_model_from_config(
            joinpath(@__DIR__, "MITgcm_2deg8_abiotic.yaml"), "abiotic_O2", modelpars=extrapars)

    elseif baseconfig == "PO4MMbase"

            model = PB.create_model_from_config(
                joinpath(@__DIR__, "MITgcm_2deg8_COPDOM.yaml"), "PO4MMbase", modelpars=extrapars)
    
    elseif baseconfig == "PO4MMbaseECCO"

            model = PB.create_model_from_config(
                joinpath(@__DIR__, "MITgcm_ECCO_COPDOM.yaml"), "PO4MMbase", modelpars=extrapars)

    elseif baseconfig == "PO4MMcarbSCH4"
            model = PB.create_model_from_config(
                    joinpath(@__DIR__, "MITgcm_2deg8_COPDOM.yaml"), "PO4MMcarbSCH4", modelpars=extrapars)
        
    else
        error("unrecognized baseconfig='$(baseconfig)'")
    end

    return model
end




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

