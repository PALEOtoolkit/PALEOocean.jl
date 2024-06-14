using Logging
using DiffEqBase
using Sundials

using Plots

import PALEOboxes as PB
import PALEOmodel
import PALEOocean
import PALEOsediment

global_logger(ConsoleLogger(stderr,Logging.Info))


include("config_ocean_shelf1D_expts.jl")
include("plot_shelf.jl")
include("../atmreservoirreaction.jl") # ReactionReservoirAtm

transport_dir = "S2P3_transport_20240614" # folder containing S2P3 physical variables output collated to netcdf files

################
# configure model options and spinup
####################

# population-based phytoplankton model + sediment
model = PB.create_model_from_config(
    joinpath(@__DIR__, "PALEO_examples_shelf1Dsed_cfg.yaml"),
    "shelf1D_P_O2_S_CH4_sed";
    modelpars=Dict(
        "phys_file"=>joinpath(transport_dir, "S2P3_depth80_m2amp04_phys.nc"),
        "surf_file"=>joinpath(transport_dir, "S2P3_depth80_m2amp04_surf.nc"),
    )
)

config_ocean_shelf1D_expts(model, ["baseline"]);

do_test = true
do_spinup = false
do_pickup = false

# do_test = false
# do_spinup = true
# do_pickup = true

# output_filename = "P_O2_S_CH4_sed_10yr_fastH2SO2_20210303"
# output_filename = "P_O2_S_CH4_sed_50yr_20210304"
tspinup=10.0
output_filename = "P_O2_S_CH4_sed_10yr_20240614"

if do_test 
    # short run for testing
    initial_state, modeldata = PALEOmodel.initialize!(model)

    # Check initial derivative:
    initial_deriv = similar(initial_state)
    PALEOmodel.SolverFunctions.ModelODE(modeldata)(initial_deriv, initial_state, nothing, 0.0)
    # Check Jacobian:
    jac, jac_prototype = PALEOmodel.JacobianAD.jac_config_ode(:ForwardDiffSparse, model, initial_state, modeldata, 0.0)
    Jtest = copy(jac_prototype)
    jac(Jtest, initial_state, nothing, 0.0)

    tspan=(0, 1.0) 

    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    sol = PALEOmodel.ODE.integrateForwardDiff(
        paleorun, initial_state, modeldata, tspan, 
        solvekwargs=(
            reltol=1e-3, 
            abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all),
            maxiters=1000000000, 
            saveat=collect(0:0.01:10.0),
        )
    ) 

    colT=collect(range(tspan[1], stop=tspan[end], step=0.5))
 
end

if do_spinup
    initial_state, modeldata = PALEOmodel.initialize!(model)

    tspan=(0, tspinup) 
    saveat=[[0.0, 0.1, 1.0]; collect(range(10.0, stop=tspinup, step=10.0))]

    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    sol = PALEOmodel.ODE.integrateForwardDiff(
        paleorun, initial_state, modeldata, tspan, 
        solvekwargs=(
            reltol=1e-5, 
            abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all),
            maxiters=1000000000,
            saveat=saveat
        )
    ) 

    # sol = PALEOmodel.ODE.integrateDAEForwardDiff(paleorun, initial_state, modeldata, tspan, solvekwargs=(saveat=1/8760, reltol=1e-5, maxiters=200000))  # first run includes JIT time

    # sol = PALEOmodel.ODE.integrate(paleorun, initial_state, modeldata, tspan, solvekwargs=(reltol=1e-5,))  # first run includes JIT time

    isempty(output_filename) || PALEOmodel.OutputWriters.save_netcdf(paleorun.output, output_filename)

    colT = saveat
end

if do_pickup
    pickup_output = PALEOmodel.OutputWriters.load_netcdf!(PALEOmodel.OutputWriters.OutputMemory(), output_filename)
    tstart = PB.get_data(pickup_output, "ocean.tmodel")[end]
    initial_state, modeldata = PALEOmodel.initialize!(model, pickup_output=pickup_output)  
   
    tspan=(tstart, tstart + 1.0) 

    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    sol = PALEOmodel.ODE.integrateForwardDiff(
        paleorun, initial_state, modeldata, tspan, 
        solvekwargs=(
            reltol=1e-5,
            abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all),
            maxiters=1000000,
            saveat=1/8760
        )
    ) 
      
    colT=collect(range(tstart, step=0.1, length=11))
end

##############################
# Plot output
###############################

# individual plots
# plotlyjs(size=(750, 565))
# pager = PALEOmodel.DefaultPlotPager()

# assemble plots onto screens with 6 subplots
gr(size=(1600, 900))
pager=PALEOmodel.PlotPager((2, 3), (legend_background_color=nothing, margin=(5, :mm)))

plot_shelf_phys(paleorun.output; pager=pager)
plot_totals_S(paleorun.output; pager=pager)
pager(:newpage)
plot_tracers_conc(
    paleorun.output; 
    tracers=["O2", "SO4", "H2S", "CH4"],
    colT=colT,
    plot_totals=true,
    pager=pager
)
plot_tracers_conc(
    paleorun.output;
    domain="sediment",
    tracers=["Corg", "O2", "SO4", "H2S", "CH4"],
    colT=colT,
    plot_totals=true,
    pager=pager
)

plot_airsea(paleorun.output; tracers=["O2"], pager=pager)
plot_biota(paleorun.output, colT=colT, pager=pager)
plot_oceanfloor(paleorun.output; pager=pager)
pager(:newpage) # flush output

