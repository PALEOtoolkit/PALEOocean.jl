using Logging
import LoggingExtras
using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOocean

global_logger(ConsoleLogger(stderr,Logging.Info))
include("config_mitgcm_expts.jl")
include("plot_mitgcm.jl")

use_threads = true

outfile_root = "MITgcm_PO4MMbaseECCO_100yr_20241220"

logfile_name = outfile_root*"_1_log.txt"
println("writing log output to ", logfile_name);  flush(stdout)
logfile = open(logfile_name, "w")
global_logger(LoggingExtras.MinLevelLogger(LoggingExtras.FileLogger(logfile, always_flush=true), Logging.Info))

model = PB.create_model_from_config(
    joinpath(@__DIR__, "MITgcm_ECCO_PO4MMbase.yaml"), "PO4MMbase";
    modelpars=Dict(
        "threadsafe"=>use_threads,
        "transport_pack_chunk_width" => 4, # use SIMD optimization for transport matrix
    ),
)
toutputs_relative = [0.0, 1.0, 10.0, 100.0]
# toutputs_relative = [0.0, 1.0]


tstep_explicit_s::Int = 43200  # s timestep to use for explicit transport
tstep_implicit_s::Int = tstep_explicit_s
tstep_explicit_yr = tstep_explicit_s/PB.Constants.k_secpyr # yr
@info "using timesteps tstep_implicit_s $tstep_implicit_s tstep_explicit_s $tstep_explicit_s tstep_explicit_yr $tstep_explicit_yr"
transportMITgcm = PB.get_reaction(model, "ocean", "transportMITgcm")
PB.setvalue!(transportMITgcm.pars.Aimp_deltat, tstep_implicit_s)
transportMITgcm = nothing  # holds large transport matrix arrays !


num_segments = 20
# num_segments = 2

toutputs=[]

if use_threads
    method_barrier = PB.reaction_method_thread_barrier(
        PALEOmodel.ThreadBarriers.ThreadBarrierAtomic("the barrier"),
        PALEOmodel.ThreadBarriers.wait_barrier
    )
else
    method_barrier = nothing
end

initial_state, modeldata = PALEOmodel.initialize!(model; method_barrier) 

for iseg in 1:num_segments
# for iseg in 20:num_segments
# iseg = 1
# if true
    @info """

    ==============================================================================
    start iseg $iseg 
    ===============================================================================
    """

    if iseg > 1
        pickup_filename = build_outfilename(outfile_root, iseg-1)
        pickup_output = PALEOmodel.OutputWriters.load_netcdf!(PALEOmodel.OutputWriters.OutputMemory(), pickup_filename)
                
        PALEOmodel.set_statevar_from_output!(modeldata, pickup_output)
        global initial_state = PALEOmodel.get_statevar(modeldata.solver_view_all)
        tstart = PB.get_data(pickup_output, "ocean.tmodel")[end]
        pickup_output = nothing        
    else
        tstart = 0.0
    end

    toutputs = toutputs_relative .+ tstart
    @info "toutputs: $toutputs"

    global paleorun = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

    if !use_threads
        @time PALEOmodel.ODEfixed.integrateEuler(paleorun, initial_state, modeldata, toutputs , tstep_explicit_yr)
    else
        cellranges = PB.Grids.get_tiled_cellranges(paleorun.model, Threads.nthreads(), "ocean")

        @time PALEOmodel.ODEfixed.integrateEulerthreads(paleorun, initial_state, modeldata, cellranges, toutputs , tstep_explicit_yr)
    end

    output_filename = build_outfilename(outfile_root, iseg)
    @info """
        ==============================================================================
        end iseg $iseg
        toutputs $toutputs 
        output_filename $output_filename
        ===============================================================================

        """
    PALEOmodel.OutputWriters.save_netcdf(paleorun.output, output_filename)
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
pager = PALEOmodel.PlotPager((2,2), (legend_background_color=nothing, margin=(5, :mm)))

plot_forcings(paleorun.output, pager=pager, lon=200.0)
pager(:newpage)
plot_abiotic_O2(paleorun.output, toutputs=toutputs, pager=pager, lon1=200.0, lon2=340.0)
pager(:newpage)
plot_PO4MMbase(
    paleorun.output,
    toutputs=toutputs,
    tbioprod=[0.5, 1.0],
    pager=pager,
)
pager(:newpage)


if !isnothing(logfile)
    close(logfile); logfile = nothing
    global_logger(ConsoleLogger(stderr, Logging.Info))
end
