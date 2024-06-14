using Logging
using DiffEqBase
using Sundials

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOocean


global_logger(ConsoleLogger(stderr,Logging.Info))

include("config_ocean_shelf1D_expts.jl")
include("plot_shelf.jl")
include("../atmreservoirreaction.jl") # ReactionReservoirAtm

transport_dir = "S2P3_transport_20240614" # folder containing S2P3 physical variables output collated to netcdf files

# population-based phytoplankton model, with sulphur and carbonate system
# NB: requires DAE solver
model = PB.create_model_from_config(
    joinpath(@__DIR__, "PALEO_examples_shelf1D_cfg.yaml"),
    "shelf1D_P_O2_S_CH4_carb";
    modelpars=Dict(
        "phys_file"=>joinpath(transport_dir, "S2P3_depth80_m2amp04_phys.nc"),
        "surf_file"=>joinpath(transport_dir, "S2P3_depth80_m2amp04_surf.nc"),
    )
)

expts = ["baseline"]; tspan=(0, 5.0)
# expts = [("setO2", 0.1)]; tspan=(0,5.0)
# expts = [("setO2", 0.0), "nooxphot", ("setH2S", 0.0 ), ("setCH4", 0.2 )]; tspan=(0,5.0)
# expts = [("setO2", 0.0), "nooxphot", ("setH2S", 0.0 ), ("setCH4", 0.1 )]; tspan=(0,5.0)
# expts = [("setO2", 0.0), ("setH2S", 0.0 ), ("setCH4", 0.0 )]; tspan=(0,5.0)
# expts = [("setO2", 0.1), ("setH2S", 0.0 ), ("setCH4", 0.0 )]; tspan=(0,5.0)
# expts = [("setO2", 0.01), ("setH2S", 0.0 ), ("setCH4", 0.0 )]; tspan=(0,5.0)

config_ocean_shelf1D_expts(model, expts)

initial_state, modeldata = PALEOmodel.initialize!(model)

# Check initial derivative:
# initial_deriv = similar(initial_state)
# PALEOmodel.SolverFunctions.ModelODE(modeldata)(initial_deriv, initial_state , nothing, 0.0)
# Check Jacobian:
# jac, jac_prototype = PALEOmodel.JacobianAD.jac_config_ode(:ForwardDiffSparse, model, initial_state, modeldata, 0.0)
# J = copy(jac_prototype)
# jac(jac_prototype, initial_state, nothing, 0.0)

paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

# first run includes JIT time
# sol = PALEOmodel.ODE.integrateForwardDiff(
#     run, initial_state, modeldata, tspan,
#     solvekwargs=(saveat=1/8760, reltol=1e-5, maxiters=1000000),
# )  

sol = PALEOmodel.ODE.integrateDAEForwardDiff(
    paleorun, initial_state, modeldata, tspan,
    solvekwargs=(saveat=1/8760, reltol=1e-5, maxiters=200000),
)

########################################
# Plot output
########################################

colT=collect(range(tspan[1], stop=tspan[end], step=0.5))

# individual plots
# plotlyjs(size=(750, 565))
# pager = PALEOmodel.DefaultPlotPager()

# assemble plots onto screens with 6 subplots
gr(size=(1600, 900))
pager=PALEOmodel.PlotPager((2, 3), (legend_background_color=nothing, ))

plot_shelf_phys(paleorun.output; pager=pager)
plot_totals_S(paleorun.output; pager=pager)
pager(:newpage)
plot_tracers_conc(
    paleorun.output; 
    tracers=["O2", "SO4", "H2S", "CH4", "TAlk", "DIC"],
    colT=colT,
    plot_totals=true,
    pager=pager
)
plot_tracers_other(
    paleorun.output; 
    tracers=["pHtot", "carbchem/OmegaAR"],
    colT=colT,
    pager=pager
)
plot_airsea(paleorun.output; tracers=["O2"], pager=pager)
plot_biota_S(paleorun.output, colT=colT, pager=pager)
pager(:newpage) # flush output




