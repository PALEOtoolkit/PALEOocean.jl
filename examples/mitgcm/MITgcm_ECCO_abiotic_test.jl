using Logging

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOocean

global_logger(ConsoleLogger(stderr,Logging.Info))
include("config_mitgcm_expts.jl")
include("plot_mitgcm.jl")

use_threads = true
do_benchmarks = false

model = PB.create_model_from_config(
    joinpath(@__DIR__, "MITgcm_ECCO_COPDOM.yaml"), "PO4MMbase"
)

toutputs = [0.0, 1.0] # , 10.0, 100.0]

transportMITgcm = PB.get_reaction(model, "ocean", "transportMITgcm")
tstep = transportMITgcm.par_Aimp_deltat[]/PB.Constants.k_secpyr
transportMITgcm = nothing

@info "using tstep=$tstep yr"

output_filename = ""

pickup_output = nothing
initial_state, modeldata = PALEOmodel.initialize!(model, threadsafe=use_threads, pickup_output=pickup_output)  
pickup_output = nothing

paleorun = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

if !use_threads
    @time PALEOmodel.ODEfixed.integrateEuler(paleorun, initial_state, modeldata, toutputs , tstep)
else
    # if Threads.nthreads() == 4 
    #     # indices are slightly uneven to equalize the number of active (mostly ocean) cells per thread
    #     tiles = [(1:124, :, :), (125:202, :, :), (203:275, :, :), (276:360, :, :)]  # 4 threads
    # elseif Threads.nthreads() == 8
    #     # indices are slightly uneven to equalize the number of active (mostly ocean) cells per thread
    #     tiles = [(1:62, :, :), (63:120, :, :), 
    #             (121:167, :, :),  (168:200, :, :),
    #             (201:232, :, :), (233:275, :, :),
    #             (276:325, :, :) , (326:360, :, :)]  # 8 threads
    # elseif Threads.nthreads() == 16
    #     # indices are slightly uneven to equalize the number of active (mostly ocean) cells per thread
    #     tiles = [(1:31, :, :), (32:62, :, :), (63:88, :, :),   (89:120, :, :), 
    #             (121:148, :, :), (149:167, :, :), (168:184, :, :),  (185:199, :, :),
    #             (200:215, :, :), (216:230, :, :), (231:250, :, :), (251:274, :, :),
    #             (275:304, :, :), (305:325, :, :),  (326:340, :, :),  (341:360, :, :),]  # 16 threads

    # else
    #     error("no config for nthreads=$(Threads.nthreads())")
    # end

    # cellranges = PB.Grids.get_tiled_cellranges(paleorun.model, tiles)  # vector of vectors, 1 per tile

    cellranges = PB.Grids.get_tiled_cellranges(paleorun.model, Threads.nthreads(), "ocean")

    for (t,tc) in enumerate(cellranges)
        numcells = 0
        for c in tc
            numcells += length(c.indices)
        end
        println("thread $t numcells $numcells")
    end
    @time PALEOmodel.ODEfixed.integrateEulerthreads(paleorun, initial_state, modeldata, cellranges, toutputs , tstep)
end

isempty(output_filename) || PALEOmodel.OutputWriters.save_jld2(paleorun.output, output_filename)


show(PB.show_variables(paleorun.model), allrows=true)
println()



