using Logging

using Plots; plotlyjs(size=(750, 565))

import PALEOboxes as PB
import PALEOmodel
import PALEOocean

global_logger(ConsoleLogger(stderr,Logging.Info))
include("config_mitgcm_expts.jl")

use_threads = true

# model = PB.create_model_from_config(
#     joinpath(@__DIR__, "MITgcm_2deg8_abiotic.yaml"), "abiotic_O2"
# )
# toutputs = [0, 0.25, 0.5, 0.75, 1.0] #, 10.0]


# model = PB.create_model_from_config(
#     joinpath(@__DIR__, "MITgcm_2deg8_COPDOM.yaml"), "PO4MMbase"
# )
# toutputs = [0, 0.25, 0.5, 0.75, 1.0] #, 10.0, 99.5, 100.0] #, 1000.0, 1000.5]


model = PB.create_model_from_config(
    joinpath(@__DIR__, "MITgcm_2deg8_COPDOM.yaml"), "PO4MMcarbSCH4";
    modelpars=Dict("threadsafe"=>use_threads),
)
toutputs = [0, 1.0] # , 10.0, 100.0, 1000.0, 1999.5, 2000.0, 2999.5, 3000.0]
# start with low oxygen to test marine sulphur system
PB.set_variable_attribute!(model, "atm", "O2", :initial_value, 0.1*3.71e19)
PB.set_variable_attribute!(model, "ocean", "O2", :initial_value, 0.1*0.2054)


transportMITgcm = PB.get_reaction(model, "ocean", "transportMITgcm")
tstep = transportMITgcm.par_Aimp_deltat[]/PB.Constants.k_secpyr

@info "using tstep=$tstep yr"

if use_threads
    method_barrier = PB.reaction_method_thread_barrier(
        PALEOmodel.ThreadBarriers.ThreadBarrierAtomic("the barrier"),
        PALEOmodel.ThreadBarriers.wait_barrier
    )
else
    method_barrier = nothing
end

initial_state, modeldata = PALEOmodel.initialize!(model; method_barrier)  
    
run = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

if !use_threads
    @time PALEOmodel.ODEfixed.integrateEuler(run, initial_state, modeldata, toutputs , tstep)
else
    #=
    if Threads.nthreads() == 4 
        # indices are slightly uneven to equalize the number of active (mostly ocean) cells per thread
        tiles = [(1:44, :, :), (45:72, :, :), (73:98, :, :), (99:128, :, :)]  # 4 threads
    elseif Threads.nthreads() == 8
        # indices are slightly uneven to equalize the number of active (mostly ocean) cells per thread
        tiles = [(1:23, :, :), (24:44, :, :), (45:61, :, :),  (62:72, :, :), (73:83, :, :), (84:98, :, :), (99:116, :, :) , (117:128, :, :)]  # 8 threads
    elseif Threads.nthreads() == 16
        # indices are slightly uneven to equalize the number of active (mostly ocean) cells per thread
        tiles = [(1:12, :, :), (13:22, :, :), (23:32, :, :), (33:44, :, :), 
                (45:55, :, :), (56:61, :, :),  (62:67, :, :),  (68:72, :, :),
                (73:77, :, :), (78:83, :, :), (84:90, :, :), (91:98, :, :), 
                (99:109, :, :) ,  (110:116, :, :) , (117:121, :, :), (122:128, :, :),]  # 16 threads

    else
        error("no config for nthreads=$(Threads.nthreads())")
    end

    cellranges = PB.Grids.get_tiled_cellranges(run.model, tiles)  # vector of vectors, 1 per tile
=#
    cellranges = PB.Grids.get_tiled_cellranges(run.model, Threads.nthreads(), "ocean")

    @time PALEOmodel.ODEfixed.integrateEulerthreads(run, initial_state, modeldata, cellranges, toutputs , tstep)
end


show(PB.show_variables(run.model), allrows=true)
println()

