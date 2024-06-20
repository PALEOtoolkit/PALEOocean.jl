using Logging
using DiffEqBase
using Sundials

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOocean


global_logger(ConsoleLogger(stderr,Logging.Info))

include("TransportExamples.jl")


###########################
# create and initialize
###########################

model = PB.create_model_from_config(
    joinpath(@__DIR__, "TransportExamples_cfg.yaml"), "example_diffuse"
)

initial_state, modeldata = PALEOmodel.initialize!(model)

# bodge an updated initial_state for testing
ocean_T = PB.get_data(PB.get_variable(model, "ocean.T"), modeldata) # model data array for tracer T
ocean_T[50] = 2e14*20.0*1.0  # ~ 1 mol m-3 in middle cell of first column
ocean_T[100+5] = 2e13*200.0*1.0  # ~ 1 mol m-3 in middle cell of second column
initial_state = PALEOmodel.get_statevar(modeldata.solver_view_all)

##############################
# integrate as an ODE
############################
tspan=(0.0, 5e3)

paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

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

# total
pager(plot(title="Total T", paleorun.output, ["ocean.T_total"]; ylabel="T (mol)",))

# line plots at specified times
columns = [:a, :b]
tcol = collect(0.0:100.0:1000)  # times at which to plot columns
for col in columns
    pager(
        plot(title="Ocean T_conc column :$col", paleorun.output, "ocean.T_conc", ( tmodel=tcol, column=col);
            swap_xy=true, xaxis=:log10, xlim=10 .^(-2.5, 0), labelattribute=:filter_records)
    )
end
pager(:newpage)

# heatmaps vs time
pager=PALEOmodel.PlotPager((2, 1), (legend_background_color=nothing, margin=(5, :mm)))
for col in columns
    d = PALEOmodel.get_array(paleorun.output, "ocean.T_conc"; column=col)
    d.values .= log10.(d.values)
    pager(
        heatmap(title="Ocean T_conc column :$col", d; clims=(-2.5, 0), xlim=(0, 1000)) # set scale for visibility
    )
end

pager(:newpage) # flush output

