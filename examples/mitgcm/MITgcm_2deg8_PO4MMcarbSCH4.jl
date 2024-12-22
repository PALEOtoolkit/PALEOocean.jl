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
use_split   = false

# transport matrices have deltaT=1200 s (timestep used to run the ocean model)
tstep_explicit_s::Int = 86400 # s timestep to use for explicit transport
n_inner = 8         # if use_split, number of implicit timesteps per explicit timestep
# n_inner = 72        # max value, gives implicit timestep 1200 s

pickup_output = nothing

# output_filename = "MITgcm_PO4MMcarbSCH42deg8_2000yr_20241222"
output_filename = ""

if !isempty(output_filename)
    logfile_name = output_filename*"_log.txt"
    println("writing log output to ", logfile_name);  flush(stdout)
    logfile = open(logfile_name, "w")
    global_logger(LoggingExtras.MinLevelLogger(LoggingExtras.FileLogger(logfile, always_flush=true), Logging.Info))
else
    logfile = nothing
end

model = PB.create_model_from_config(
    joinpath(@__DIR__, "MITgcm_2deg8_PO4MMcarbSCH4.yaml"), "PO4MMcarbSCH4";
    modelpars=Dict(
        "threadsafe"=>use_threads,
        "carbchem_simd_width"=>"FP32P8", # use SIMD optimisation for carbonate chemistry
        "transport_pack_chunk_width" => 4, # use SIMD optimization for transport matrix
        # "transport_pack_chunk_width" => 8, # use SIMD optimization for transport matrix
    ),
)

toutputs = [0, 1.0, 10.0] # short run for testing
# toutputs = [0, 1.0, 10.0, 100.0, 1000.0, 1999.5, 2000.0,]

# start with low oxygen to test marine sulphur system
PB.set_variable_attribute!(model, "atm", "O2", :initial_value, 0.1*3.71e19)
PB.set_variable_attribute!(model, "ocean", "O2", :initial_value, 0.1*0.2054)

# configure timestepping
tstep_explicit_yr = tstep_explicit_s/PB.Constants.k_secpyr # yr
if use_split
    tstep_implicit_s::Int = tstep_explicit_s / n_inner
else
    tstep_implicit_s::Int = tstep_explicit_s 
end
@info "using timesteps tstep_implicit_s $tstep_implicit_s tstep_explicit_s $tstep_explicit_s tstep_explicit_yr $tstep_explicit_yr"
transportMITgcm = PB.get_reaction(model, "ocean", "transportMITgcm")
PB.setvalue!(transportMITgcm.pars.Aimp_deltat, tstep_implicit_s)
transportMITgcm = nothing # holds large transport matrix arrays !


if use_threads
    method_barrier = PB.reaction_method_thread_barrier(
        PALEOmodel.ThreadBarriers.ThreadBarrierAtomic("the barrier"),
        PALEOmodel.ThreadBarriers.wait_barrier;
        operatorID = [1], # if use_split = true, only operatorID 1 (explict transport matrix) has dependency between tiles
    )
else
    method_barrier = nothing
end


initial_state, modeldata = PALEOmodel.initialize!(model; method_barrier, pickup_output)  
pickup_output = nothing

paleorun = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

if !use_threads
    if use_split       
        @info "using tstep_outer=$tstep_explicit_yr (yr), n_inner $n_inner, tstep_inner $(tstep_explicit_yr/n_inner) yr"
        cellrange_outer = PB.create_default_cellrange(paleorun.model, operatorID=1)
        cellrange_inner = PB.create_default_cellrange(paleorun.model, operatorID=2)
        
        @time PALEOmodel.ODEfixed.integrateSplitEuler(
            paleorun, initial_state, modeldata, toutputs, tstep_explicit_yr, n_inner,
            cellranges_outer=cellrange_outer,
            cellranges_inner=cellrange_inner
        )
    else
        @info "using tstep=$tstep_explicit_yr (yr)"
        @time PALEOmodel.ODEfixed.integrateEuler(paleorun, initial_state, modeldata, toutputs, tstep_explicit_yr)
    end
else
    if use_split
        @info "using tstep_outer=$tstep_explicit_yr (yr), n_inner $n_inner, tstep_inner $(tstep_explicit_yr/n_inner) yr"
        cellranges_outer = PB.Grids.get_tiled_cellranges(paleorun.model, Threads.nthreads(), "ocean", operatorID=1)
        cellranges_inner = PB.Grids.get_tiled_cellranges(paleorun.model, Threads.nthreads(), "ocean", operatorID=2)
        @time PALEOmodel.ODEfixed.integrateSplitEulerthreads(
            paleorun, initial_state, modeldata, toutputs , tstep_explicit_yr, n_inner;
            cellranges_outer=cellranges_outer, cellranges_inner=cellranges_inner
        )
    else
        @info "using tstep=$tstep_explicit_yr (yr)"
        cellranges = PB.Grids.get_tiled_cellranges(paleorun.model, Threads.nthreads(), "ocean")
        @time PALEOmodel.ODEfixed.integrateEulerthreads(paleorun, initial_state, modeldata, cellranges, toutputs , tstep_explicit_yr)
    end
end

isempty(output_filename) || PALEOmodel.OutputWriters.save_netcdf(paleorun.output, output_filename)

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



if !isnothing(logfile)
    close(logfile); logfile = nothing
    global_logger(ConsoleLogger(stderr, Logging.Info))
end