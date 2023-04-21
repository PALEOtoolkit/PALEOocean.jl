using Logging

using Plots; plotlyjs(size=(750, 565))

import PALEOboxes as PB
import PALEOmodel
import PALEOocean

global_logger(ConsoleLogger(stderr,Logging.Info))
include("config_mitgcm_expts.jl")
include("plot_mitgcm.jl")

use_threads = false
use_split   = false
n_inner = 2

# model = config_mitgcm_expts("PO4MMbase", ""); toutputs = [0.0, 1.0, 10.0, 100.0, 995,0, 1000.0, 1999.5, 2000.0, 2999.5, 3000.0] #, 1000.0, 1000.5]

model = PB.create_model_from_config(
    joinpath(@__DIR__, "MITgcm_2deg8_COPDOM.yaml"), "PO4MMbase"
)

toutputs = [0.0, 0.25, 0.5, 0.75, 1.0, 10.0] 
# toutputs = [0.0, 1.0, 10.0, 100.0, 995,0, 1000.0, 1999.5, 2000.0, 2999.5, 3000.0] #, 1000.0, 1000.5]

transportMITgcm = PB.get_reaction(model, "ocean", "transportMITgcm")
tstep_imp = transportMITgcm.pars.Aimp_deltat[]/PB.Constants.k_secpyr


output_filename = ""
# output_filename = "MITgcm_PO4MMbase2deg8_3000yr_20210202"

pickup_output = nothing
initial_state, modeldata = PALEOmodel.initialize!(model, threadsafe=use_threads, pickup_output=pickup_output)  
pickup_output = nothing

paleorun = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

if !use_threads
    if use_split       
        @info "using tstep_outer=$n_inner x $tstep_imp yr"
        cellrange_outer = PB.create_default_cellrange(paleorun.model, operatorID=1)
        cellrange_inner = PB.create_default_cellrange(paleorun.model, operatorID=2)
        
        @time PALEOmodel.ODEfixed.integrateSplitEuler(
            paleorun, initial_state, modeldata, toutputs, tstep_imp*n_inner, n_inner,
            cellrange_outer=cellrange_outer,
            cellrange_inner=cellrange_inner
        )
    else
        @info "using tstep=$tstep_imp yr"
        @time PALEOmodel.ODEfixed.integrateEuler(paleorun, initial_state, modeldata, toutputs, tstep_imp)
    end
else
    # Threads.nthreads() == 4 || error("use_threads requires 4 threads, Threads.nthreads()=", Threads.nthreads())

    # indices are slightly uneven to equalize the number of active (mostly ocean) cells per thread
    # tiles = [(1:44, :, :), (45:72, :, :), (73:98, :, :), (99:128, :, :)]  # 4 threads

    # cellranges = PB.Grids.get_tiled_cellranges(paleorun.model, tiles)  # vector of vectors, 1 per tile

    cellranges = PB.Grids.get_tiled_cellranges(paleorun.model, Threads.nthreads(), "ocean")

    if use_split
        @info "using tstep_outer=$n_inner x $tstep_imp yr"
        cellranges_outer = PB.Grids.get_tiled_cellranges(paleorun.model, Threads.nthreads(), "ocean", operatorID=1)
        cellranges_inner = PB.Grids.get_tiled_cellranges(paleorun.model, Threads.nthreads(), "ocean", operatorID=2)
        @time PALEOmodel.ODEfixed.integrateSplitEulerthreads(paleorun, initial_state, modeldata, toutputs , tstep_imp*n_inner, n_inner,
                cellranges_outer=cellranges_outer, cellranges_inner=cellranges_inner)
    else
        @info "using tstep=$tstep_imp yr"
        @time PALEOmodel.ODEfixed.integrateEulerthreads(paleorun, initial_state, modeldata, cellranges, toutputs , tstep_imp)
    end
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
pager = PALEOmodel.PlotPager((2,2), (legend_background_color=nothing, ))

plot_forcings(paleorun.output, pager=pager)
pager(:newpage)
plot_abiotic_O2(paleorun.output, toutputs=toutputs, pager=pager)
pager(:newpage)
plot_PO4MMbase(
    paleorun.output,
    toutputs=toutputs,
    tbioprod=[0.5, 1.0],
    pager=pager,
)
pager(:newpage)
