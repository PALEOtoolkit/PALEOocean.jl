
using Logging

import PALEOboxes as PB
import PALEOmodel
import PALEOocean


global_logger(ConsoleLogger(stderr, Logging.Info))

# Skeleton ocean-atmosphere configuration, with no biogeochemistry

model = PB.create_model_from_config(
   joinpath(@__DIR__, "PALEO_examples_oceanskeleton_cfg.yaml"), "ocean_skeleton_COP")

initial_state, modeldata = PALEOmodel.initialize!(model)

# call ODE function to check derivative
println("initial_state", initial_state)
initial_deriv = similar(initial_state)
PALEOmodel.SolverFunctions.ModelODE(modeldata)(initial_deriv, initial_state , nothing, 0.0)
println("initial_deriv", initial_deriv)

println()
println("show_variables:")
show(PB.show_variables(model); allrows=true)
