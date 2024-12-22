using Logging

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

model = PB.create_model_from_config(
    joinpath(@__DIR__, "MITgcm_2deg8_PO4MMbase.yaml"), "PO4MMbase";
    modelpars=Dict(
        "threadsafe"=>use_threads,
        "transport_pack_chunk_width" => 4, # use SIMD optimization for transport matrix
    ),
)

toutputs = [0.0, 0.25, 0.5, 0.75, 1.0, 10.0] 
# toutputs = [0.0, 1.0, 10.0, 100.0, 995,0, 1000.0, 1999.5, 2000.0, 2999.5, 3000.0] #, 1000.0, 1000.5]

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


output_filename = ""
# output_filename = "MITgcm_PO4MMbase2deg8_3000yr_20210202"

if use_threads
    method_barrier = PB.reaction_method_thread_barrier(
        PALEOmodel.ThreadBarriers.ThreadBarrierAtomic("the barrier"),
        PALEOmodel.ThreadBarriers.wait_barrier;
        operatorID = [1], # if use_split = true, only operatorID 1 (explict transport matrix) has dependency between tiles
    )
else
    method_barrier = nothing
end

pickup_output = nothing
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
    cellranges = PB.Grids.get_tiled_cellranges(paleorun.model, Threads.nthreads(), "ocean")

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
plot_abiotic_O2(paleorun.output, toutputs=toutputs, pager=pager)
pager(:newpage)
plot_PO4MMbase(
    paleorun.output,
    toutputs=toutputs,
    tbioprod=[0.5, 1.0],
    pager=pager,
)
pager(:newpage)
