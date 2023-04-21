using Logging

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOocean

global_logger(ConsoleLogger(stderr,Logging.Info))
include("config_mitgcm_expts.jl")
include("plot_mitgcm.jl")

use_threads = false


model = PB.create_model_from_config(
    joinpath(@__DIR__, "MITgcm_2deg8_abiotic.yaml"), "abiotic_O2"
)

toutputs_relative = [0, 0.25, 0.5, 0.75, 1.0] #, 10.0]

transportMITgcm = PB.get_reaction(model, "ocean", "transportMITgcm")
tstep = transportMITgcm.pars.Aimp_deltat[]/PB.Constants.k_secpyr

@info "using tstep=$tstep yr"

num_segments = 1
outfile_root = ""
# num_segments = 20
# outfile_root = "MITgcm_PO4MMbase2deg8_100yr_20210201"

toutputs = []

for iseg in 1:num_segments
    if iseg > 1
        !isempty(outfile_root) || error("outfile_root is empty")
        pickup_filename = build_outfilename(outfile_root, iseg-1)
        pickup_output = PALEOmodel.OutputWriters.load_jld2!(PALEOmodel.OutputWriters.OutputMemory(), pickup_filename)
       
        tstart = PB.get_data(pickup_output, "ocean.tmodel")[end] 
    
    else
        pickup_output = nothing
        tstart = 0.0
    end

    global initial_state
    global modeldata
    global toutputs

    initial_state, modeldata = PALEOmodel.initialize!(model, threadsafe=use_threads, pickup_output=pickup_output) 
    pickup_output = nothing

    toutputs = toutputs_relative .+ tstart

    global paleorun = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

    if !use_threads
        @time PALEOmodel.ODEfixed.integrateEuler(paleorun, initial_state, modeldata, toutputs , tstep)
    else
        # Threads.nthreads() == 4 || error("use_threads requires 4 threads, Threads.nthreads()=", Threads.nthreads())

        # # indices are slightly uneven to equalize the number of active (mostly ocean) cells per thread
        # tiles = [(1:44, :, :), (45:72, :, :), (73:98, :, :), (99:128, :, :)]  # 4 threads
    
        # cellranges = PB.Grids.get_tiled_cellranges(paleorun.model, tiles)  # vector of vectors, 1 per tile

        cellranges = PB.Grids.get_tiled_cellranges(run.model, Threads.nthreads(), "ocean")

        @time PALEOmodel.ODEfixed.integrateEulerthreads(paleorun, initial_state, modeldata, cellranges, toutputs , tstep)
    end

    if !isempty(outfile_root)
        output_filename = build_outfilename(outfile_root, iseg)
        PALEOmodel.OutputWriters.save_jld2(paleorun.output, output_filename)
    end
end

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


