
using Logging
using DiffEqBase
using Sundials

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOocean
import PALEOcopse


global_logger(ConsoleLogger(stderr, Logging.Info))

include("config_ocean3box_expts.jl")
include("plot_ocean_3box.jl")


model = config_ocean3box_expts("oaopencarb", ["killbio", "lowO2"]); tspan=(-10e6, 10e6) # tspan=(-10e6,1000.0) # 

initial_state, modeldata = PALEOmodel.initialize!(model)
statevar_norm = PALEOmodel.get_statevar_norm(modeldata.solver_view_all)

# call ODE function to check derivative
println("initial_state", initial_state)
println("statevar_norm", statevar_norm)
initial_deriv = similar(initial_state)
PALEOmodel.SolverFunctions.ModelODE(modeldata)(initial_deriv, initial_state , nothing, 0.0)
println("initial_deriv", initial_deriv)

run = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

# With `killbio` H2S goes to zero, so this provides a test case for solvers `abstol` handling
# (without this option, solver will fail or take excessive steps as it attempts to solve H2S for noise) 

# Solve as DAE with sparse Jacobian
PALEOmodel.ODE.integrateDAEForwardDiff(
   run, initial_state, modeldata, tspan,
   alg=IDA(linear_solver=:KLU),
   solvekwargs=(
      abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all),
      save_start=false
   )
)

# Solve as ODE with Jacobian (OK if no carbonate chem or global temperature)
# sol = PALEOmodel.ODE.integrateForwardDiff(run, initial_state, modeldata, tspan, alg=CVODE_BDF(linear_solver=:KLU))
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

plot_totals(run.output; species=["C", "TAlk", "TAlkerror", "O2", "S", "P"], pager=pager)
plot_ocean_tracers(
    run.output; 
    tracers=["TAlk_conc", "DIC_conc", "temp", "pHtot", "O2_conc", "SO4_conc", "H2S_conc", "P_conc", 
        "SO4_delta", "H2S_delta", "pHtot", "OmegaAR"],
    pager=pager
)
plot_oaonly_abiotic(run.output; pager=pager)
plot_carb_open(run.output; pager=pager)
pager(:newpage) # flush output
