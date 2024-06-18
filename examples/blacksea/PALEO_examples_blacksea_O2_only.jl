using Logging
using DiffEqBase
using Sundials

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOocean


global_logger(ConsoleLogger(stderr,Logging.Info))

include("config_ocean_blacksea_expts.jl")
include("plot_ocean_blacksea.jl")

# abiotic atmosphere-ocean, O2 only
model = PB.create_model_from_config(
    joinpath(@__DIR__, "PALEO_examples_blacksea_cfg.yaml"), "blacksea_abiotic_O2",
    modelpars=Dict(
        "matdir"=>joinpath(@__DIR__, "../romglb/romaniello2010_transport") # assume zip file has been downloaded and unpacked in this subfolder
    )
)

config_ocean_blacksea_expts(model, ["baseline"]); tspan=(0, 1e4)

initial_state, modeldata = PALEOmodel.initialize!(model)

paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

# No carbonate system, so can integrate as an ODE
PALEOmodel.ODE.integrateForwardDiff(
    paleorun, initial_state, modeldata, tspan,
    solvekwargs=(reltol=1e-5,),
)


########################################
# Plot output
########################################

# individual plots
# plotlyjs(size=(750, 565))
# pager = PALEOmodel.DefaultPlotPager()

# assemble plots onto screens with 6 subplots
gr(size=(1200, 900))
pager=PALEOmodel.PlotPager((2, 3), (legend_background_color=nothing, margin=(5, :mm)))

plot_totals(paleorun.output; species=["T", "O2"], pager=pager)
plot_airsea(paleorun.output; species=["O2"], pager=pager)
plot_ocean_tracers(
    paleorun.output; 
    tracers=["insol", "T_conc", "O2_conc"],
    tcol=[-Inf, 10.0, 100.0, 1000.0, 1e4, 1e5],
    pager=pager
)
pager(:newpage) # flush output
