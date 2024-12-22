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
    joinpath(@__DIR__, "MITgcm_2deg8_abiotic.yaml"), "abiotic_O2";
    modelpars=Dict(
        "threadsafe"=>use_threads,
        "transport_pack_chunk_width" => 4, # use SIMD optimization for transport matrix
    ),
)

toutputs = [0, 0.25, 0.5, 0.75, 1.0] #, 10.0]

tstep_explicit_s::Int = 86400 # s timestep to use for explicit transport
tstep_implicit_s::Int = tstep_explicit_s
tstep_explicit_yr = tstep_explicit_s/PB.Constants.k_secpyr # yr
@info "using timesteps tstep_implicit_s $tstep_implicit_s tstep_explicit_s $tstep_explicit_s tstep_explicit_yr $tstep_explicit_yr"
transportMITgcm = PB.get_reaction(model, "ocean", "transportMITgcm")
PB.setvalue!(transportMITgcm.pars.Aimp_deltat, tstep_implicit_s)
transportMITgcm = nothing # holds large transport matrix arrays !

output_filename = ""

if use_threads
    method_barrier = PB.reaction_method_thread_barrier(
        PALEOmodel.ThreadBarriers.ThreadBarrierAtomic("the barrier"),
        PALEOmodel.ThreadBarriers.wait_barrier
    )
else
    method_barrier = nothing
end

initial_state, modeldata = PALEOmodel.initialize!(model; method_barrier) 

paleorun = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

if !use_threads
    @time PALEOmodel.ODEfixed.integrateEuler(paleorun, initial_state, modeldata, toutputs , tstep_explicit_yr)
else
    cellranges = PB.Grids.get_tiled_cellranges(run.model, Threads.nthreads(), "ocean")

    @time PALEOmodel.ODEfixed.integrateEulerthreads(paleorun, initial_state, modeldata, cellranges, toutputs , tstep_explicit_yr)
end

if !isempty(output_filename)
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

plot_forcings(paleorun.output, pager=pager)
pager(:newpage)
plot_abiotic_O2(paleorun.output, toutputs=toutputs, pager=pager)
pager(:newpage)


