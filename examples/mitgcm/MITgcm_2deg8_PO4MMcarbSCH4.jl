using Logging

using Plots; plotlyjs(size=(750, 565))

import PALEOboxes as PB
import PALEOmodel
import PALEOocean

global_logger(ConsoleLogger(stderr,Logging.Info))
include("config_mitgcm_expts.jl")
include("plot_mitgcm.jl")

use_threads = true

model = PB.create_model_from_config(
    joinpath(@__DIR__, "MITgcm_2deg8_COPDOM.yaml"), "PO4MMcarbSCH4"
)

toutputs = [0, 1.0, 10.0] # , 100.0, 1000.0, 1999.5, 2000.0, 2999.5, 3000.0]
# toutputs = [0, 1.0, 10.0, 100.0, 1000.0, 1999.5, 2000.0,]

# start with low oxygen to test marine sulphur system
PB.set_variable_attribute!(model, "atm", "O2", :initial_value, 0.1*3.71e19)
PB.set_variable_attribute!(model, "ocean", "O2", :initial_value, 0.1*0.2054)

transportMITgcm = PB.get_reaction(model, "ocean", "transportMITgcm")
tstep = transportMITgcm.pars.Aimp_deltat[]/PB.Constants.k_secpyr
@info "using tstep=$tstep yr"

# output_filename = "MITgcm_PO4MMcarbSCH42deg8FP64_3000yr_20210210"
# output_filename = "MITgcm_PO4MMcarbSCH42deg8FP64_2000yr_20230422"
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
    #     tiles = [(1:44, :, :), (45:72, :, :), (73:98, :, :), (99:128, :, :)]  # 4 threads
    # elseif Threads.nthreads() == 8
    #     # indices are slightly uneven to equalize the number of active (mostly ocean) cells per thread
    #     tiles = [(1:23, :, :), (24:44, :, :), (45:61, :, :),  (62:72, :, :), (73:83, :, :), (84:98, :, :), (99:116, :, :) , (117:128, :, :)]  # 8 threads
    # elseif Threads.nthreads() == 16
    #     # indices are slightly uneven to equalize the number of active (mostly ocean) cells per thread
    #     tiles = [(1:12, :, :), (13:22, :, :), (23:32, :, :), (33:44, :, :), 
    #             (45:55, :, :), (56:61, :, :),  (62:67, :, :),  (68:72, :, :),
    #             (73:77, :, :), (78:83, :, :), (84:90, :, :), (91:98, :, :), 
    #             (99:109, :, :) ,  (110:116, :, :) , (117:121, :, :), (122:128, :, :),]  # 16 threads

    # else
    #     error("no config for nthreads=$(Threads.nthreads())")
    # end

    # cellranges = PB.Grids.get_tiled_cellranges(paleorun.model, tiles)  # vector of vectors, 1 per tile

    cellranges = PB.Grids.get_tiled_cellranges(paleorun.model, Threads.nthreads(), "ocean")

    @time PALEOmodel.ODEfixed.integrateEulerthreads(paleorun, initial_state, modeldata, cellranges, toutputs , tstep)
end

isempty(output_filename) || PALEOmodel.OutputWriters.save_jld2(paleorun.output, output_filename)

show(PB.show_variables(paleorun.model), allrows=true)
println()

############################
# Plot 
############################

# single plots
# plotlyjs(size=(750, 565))
# pager = PALEOmodel.DefaultPlotPager()

# multiple plots per screen
gr(size=(1200, 900))
pager = PALEOmodel.PlotPager((2,2), (legend_background_color=nothing, margin=(5, :mm)))

plot_forcings(paleorun.output, pager=pager)
pager(:newpage)
plot_PO4MMbase(
    paleorun.output,
    toutputs=toutputs,
    tbioprod=[0.5, 1.0],
    pager=pager,
)
pager(:newpage)
plot_tracers(
    paleorun.output,
    tracers=["O2_conc", "SO4_conc", "H2S_conc", "CH4_conc", "SO4_delta", "H2S_delta", "CH4_delta", "TAlk_conc", "DIC_conc", "DIC_delta"],
    toutputs=[1e12],
    pager=pager
)
pager(:newpage)
plot_carbSCH4(paleorun.output, pager=pager)
pager(:newpage)



