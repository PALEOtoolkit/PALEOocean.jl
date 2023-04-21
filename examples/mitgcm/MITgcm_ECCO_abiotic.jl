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

model = config_mitgcm_expts("PO4MMbaseECCO", ""); toutputs_relative = [0.0, 1.0, 10.0, 100.0]

transportMITgcm = PB.get_reaction(model, "ocean", "transportMITgcm")
tstep = transportMITgcm.par_Aimp_deltat[]/PB.Constants.k_secpyr
transportMITgcm = nothing

@info "using tstep=$tstep yr"

num_segments = 20
outfile_root = "MITgcm_PO4MMbaseECCO_100yr_20210130"

toutputs=[]
for iseg in 1:num_segments
# for iseg in 13:num_segments
    if iseg > 1
        pickup_filename = build_outfilename(outfile_root, iseg-1)
        pickup_output = PALEOmodel.OutputWriters.load_jld2!(PALEOmodel.OutputWriters.OutputMemory(), pickup_filename)
        tstart = PB.get_data(pickup_output, "ocean.tmodel")[end]
    else
        pickup_output = nothing
        tstart = 0.0
    end
    initial_state, modeldata = PALEOmodel.initialize!(model, threadsafe=use_threads, pickup_output=pickup_output)  
    pickup_output = nothing

    toutputs = toutputs_relative .+ tstart

    global paleorun = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

    if !use_threads
        @time PALEOmodel.ODEfixed.integrateEuler(paleorun, initial_state, modeldata, toutputs , tstep)
    else
        # Threads.nthreads() == 4 || error("use_threads requires 4 threads, Threads.nthreads()=", Threads.nthreads())

        # # indices are slightly uneven to equalize the number of active (mostly ocean) cells per thread
        # # tiles = [(1:44, :, :), (45:72, :, :), (73:98, :, :), (99:128, :, :)]  # 4 threads
        # tiles = [(1:124, :, :), (125:202, :, :), (203:275, :, :), (276:360, :, :)]  # 4 threads
    
        # cellranges = PB.Grids.get_tiled_cellranges(paleorun.model, tiles)  # vector of vectors, 1 per tile

        cellranges = PB.Grids.get_tiled_cellranges(paleorun.model, Threads.nthreads(), "ocean")

        @time PALEOmodel.ODEfixed.integrateEulerthreads(paleorun, initial_state, modeldata, cellranges, toutputs , tstep)
    end

    output_filename = build_outfilename(outfile_root, iseg)
    PALEOmodel.OutputWriters.save_jld2(paleorun.output, output_filename)
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

plot_forcings(paleorun.output, pager=pager, lonidx=200)
pager(:newpage)
plot_abiotic_O2(paleorun.output, toutputs=toutputs, pager=pager, lonidx1=200, lonidx2=340)
pager(:newpage)
plot_PO4MMbase(
    paleorun.output,
    toutputs=toutputs,
    tbioprod=[0.5, 1.0],
    pager=pager,
)
pager(:newpage)


