
using Logging
using DiffEqBase
using Sundials

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOreactions
import PALEOcopse


global_logger(ConsoleLogger(stderr,Logging.Info))

include("config_PTB3box_expts.jl")
include("../plot_ocean_3box.jl")

# model = config_PTB3box_expts("Co2HOmLWCpp", ["baseline"]); tspan=(-260e6,-240e6) # tspan=(-10e6,10e6) 
# model = config_PTB3box_expts("Co2LOmHWC4pp", ["baseline"]); tspan=(-260e6,-240e6) # tspan=(-10e6,10e6) 

# Clarkson etal (2015) CO2LO scenario
model = config_PTB3box_expts("Co2LOmHWC4pp", ["Sw_2Ts", "Pp_PEes", "Lk_2", "Cia_s2", "Cib_1"]); tspan=(-260e6,-240e6) # tspan=(-10e6,10e6) 

# Clarkson etal (2015) CO2LO + kill marine biota at EP2
# model = config_PTB3box_expts("Co2LOmHWC4pp", ["Sw_2Ts", "Pp_PEes", "Lk_2", "Cia_s2", "Cib_1", "killbioEP2"]); tspan=(-260e6,-240e6)


initial_state, modeldata = PALEOmodel.initialize!(model)

# call ODE function to check derivative
initial_deriv = similar(initial_state)
PALEOmodel.SolverFunctions.ModelODE(modeldata)(initial_deriv, initial_state , nothing, 0.0)
println("initial_state", initial_state)
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
      reltol=1e-5,
      save_start=false, 
      dtmax=0.5e5, # , tstops=[-251.95e6]))
   )
)


########################################
# Plot output
########################################

# individual plots
# plotlyjs(size=(750, 565))
# pager = PALEOmodel.DefaultPlotPager()

# assemble plots onto screens with 6 subplots
gr(size=(1200, 900))

pager=PALEOmodel.PlotPager((2, 3), (xlim=(-252.15e6, -251.80e6), xflip=true, legend_background_color=nothing, ))

plot_totals(run.output; species=["C", "TAlk", "TAlkerror", "S", "P"], pager=pager)
plot_ocean_tracers(
    run.output; 
    tracers=["TAlk_conc", "DIC_conc", "temp", "pHtot", "O2_conc", "SO4_conc", "H2S_conc", "P_conc", 
        "SO4_delta", "H2S_delta", "OmegaAR", "DIC_delta", "pHtot", "BOH4_delta"],
    pager=pager
)
plot_oaonly_abiotic(run.output; pager=pager)
plot_PTB3box(run.output; pager=pager)
plot_carb_open(run.output; pager=pager)
pager(:newpage) # flush output
