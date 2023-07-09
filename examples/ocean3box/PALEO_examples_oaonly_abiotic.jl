
using Logging
using DiffEqBase
using Sundials

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOocean


global_logger(ConsoleLogger(stderr, Logging.Info))



include("config_ocean3box_expts.jl")
include("plot_ocean_3box.jl")

# Ocean-atmosphere only, test case cf Sarmiento & Toggweiler (2007) book, Fig 10.4, p436-7
# This tests the effect of air-sea exchange rate for an abiotic model.

# set k_piston below to show effect of default/fast/slow air-sea exchange rates
      
model = PB.create_model_from_config(
   joinpath(@__DIR__, "PALEO_examples_ocean3box_cfg.yaml"), "ocean3box_oaonly_abiotic_base")

config_ocean3box_expts(model, ["baseline"]); tspan=(0,1e4)
# config_ocean3box_expts(model, ["fastexchange"]); tspan=(0,1e4)
# config_ocean3box_expts(model, ["slowexchange"]); tspan=(0,1e4)


initial_state, modeldata = PALEOmodel.initialize!(model)
statevar_norm = PALEOmodel.get_statevar_norm(modeldata.solver_view_all)

# call ODE function to check derivative
println("initial_state", initial_state)
println("statevar_norm", statevar_norm)
initial_deriv = similar(initial_state)
PALEOmodel.SolverFunctions.ModelODE(modeldata)(initial_deriv, initial_state , nothing, 0.0)
println("initial_deriv", initial_deriv)

paleorun = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

# With `killbio` H2S goes to zero, so this provides a test case for solvers `abstol` handling
# (without this option, solver will fail or take excessive steps as it attempts to solve H2S for noise) 

# Solve as DAE with sparse Jacobian
PALEOmodel.ODE.integrateDAEForwardDiff(
   paleorun, initial_state, modeldata, tspan,
   alg=IDA(linear_solver=:KLU),
   solvekwargs=(
      abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all),
      # save_start=false, # fails with Sundials.jl 4.19.1
   )
)

# Solve as ODE with Jacobian (OK if no carbonate chem or global temperature)
# sol = PALEOmodel.ODE.integrateForwardDiff(paleorun, initial_state, modeldata, tspan, alg=CVODE_BDF(linear_solver=:KLU))
#    solvekwargs=(abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all),))


########################################
# Plot output
########################################

# individual plots
# plotlyjs(size=(750, 565))
# pager = PALEOmodel.DefaultPlotPager()

# assemble plots onto screens with 6 subplots
gr(size=(1200, 900))
pager=PALEOmodel.PlotPager((2, 3), (legend_background_color=nothing, ))

plot_totals(paleorun.output; species=["C", "TAlk", "TAlkerror"], pager=pager)
plot_ocean_tracers(paleorun.output; tracers=["TAlk_conc", "DIC_conc", "temp", "pHtot"], pager=pager)
plot_oaonly_abiotic(paleorun.output; pager=pager)
pager(:newpage) # flush output

#####################################
# Additional tests for solvers
########################################

# Solve as DAE without Jacobian
# PALEOmodel.ODE.integrateDAE(paleorun, initial_state, modeldata, tspan, alg=IDA(),
#    solvekwargs=(abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all), save_start=false))
