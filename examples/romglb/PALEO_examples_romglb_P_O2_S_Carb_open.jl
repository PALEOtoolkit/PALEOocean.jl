using Logging
using DiffEqBase
using Sundials

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOocean
import PALEOcopse

global_logger(ConsoleLogger(stderr,Logging.Info))

include("../atmreservoirreaction.jl")
include("SedimentationRate_dev.jl")
include("config_ocean_romglb_expts.jl")
include("plot_ocean_romglb.jl")

# Atmosphere/ocean, Ccarb, Corg, S and P burial
# Test two options for solving carbonate chemistry
model = PB.create_model_from_config(
    joinpath(@__DIR__, "PALEO_examples_romglb_cfg.yaml"), "romglb_P_O2_S_Carb_open";
    modelpars=Dict(
        "matdir"=>joinpath(@__DIR__, "romaniello2010_transport"), # assume zip file has been downloaded and unpacked in this subfolder
        # Option 1: add pHfree, TAlk to state variables and apply constraint on TAlk      
        "TAlkStateExplicit"=>true,
        # Option 2: TAlk is an implicit variable, pHfree is a state variable
        # "TAlkStateExplicit"=>false,
    ),
)

# additional configuration to set shelf and deep carbonate burial
config_ocean_romglb_expts(
    model, 
    [
        "set_carbonate_config_79box",  # configure shelf and deep carbonate burial for modern Earth
    ],
)
tspan=(0,1e5) # constraint

initial_state, modeldata = PALEOmodel.initialize!(model)

paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

# DAE solver required when using carbonate system
PALEOmodel.ODE.integrateDAEForwardDiff(
    paleorun, initial_state, modeldata, tspan,
    solvekwargs=(
        reltol=1e-5,
        dtmax=1e5,
    ),
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

plot_totals(paleorun.output; species=["P", "S"], pager=pager)
plot_airsea(paleorun.output; pager=pager)
pager(:newpage)
plot_ocean_tracers(
    paleorun.output; 
    tracers=["insol", "O2_conc", "P_conc", "H2S_conc", "DIC_conc", "TAlk_conc", "pHtot", "OmegaAR"],
    tcol=[-Inf, 10.0, 100.0, 1000.0, 1e4, 1e5],
    pager=pager
)
plot_ocean_burial(paleorun.output; pager=pager)
pager(:newpage) # flush output
