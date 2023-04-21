using Test
using BenchmarkTools

import PALEOboxes as PB

import PALEOocean
import PALEOmodel

using Plots

@testset "MITgcm" begin

skipped_testsets = [
    # "forcing", 
    # "transport"
]

do_plots = false

!("forcing" in skipped_testsets) && @testset "forcing" begin

    model = PB.create_model_from_config(
        joinpath(@__DIR__, "configoceanTMM.yaml"), "test_trspt_read")

    initial_state, modeldata = PALEOmodel.initialize!(model)

    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    # integrate to approx steady state
    @time PALEOmodel.ODEfixed.integrateEuler(paleorun, initial_state, modeldata, (0, 0.25, 0.5, 0.75, 1.0), 1/365)
    
    show(PB.show_variables(model), allrows=true)
    println()

    if do_plots
        pager = PALEOmodel.PlotPager((2, 2))
    else
        pager = PALEOmodel.NoPlotPager()
    end
   
    for t in [0.0, 0.5]
        pager(
            heatmap(paleorun.output, "oceansurface.surface_downwelling_photosynthetic_radiative_flux", (tmodel=t,), swap_xy=true),
            heatmap(paleorun.output, "oceansurface.open_area_fraction", (tmodel=t,), swap_xy=true),
        )     
    end
       
end

!("transport" in skipped_testsets) && @testset "transport" begin

    model = PB.create_model_from_config(
        joinpath(@__DIR__, "configoceanTMM.yaml"), "test_trspt_read")

    initial_state, modeldata = PALEOmodel.initialize!(model)

    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    # bodge a test for ocean tracer with single non-zero cell
    ocean_domain = PB.get_domain(model, "ocean")
    ocean_Tracer_data = PB.get_data(PB.get_variable(ocean_domain,"Tracer"), modeldata)
    Tracer_total = sum(ocean_Tracer_data)
    ocean_Tracer_data .= 0.0
    ocean_Tracer_data[1] = Tracer_total
    # update initial_state with our bodged values
    initial_state = PALEOmodel.get_statevar(modeldata.solver_view_all)

    # integrate to approx steady state
    # toutputs = [0, 0.25, 0.5, 0.75, 1.0, 10.0, 100.0, 1000.0]
    toutputs = [0, 0.25, 0.5, 0.75, 1.0]
    @time PALEOmodel.ODEfixed.integrateEuler(paleorun, initial_state, modeldata, toutputs , 1/365)
    
    show(PB.show_variables(model), allrows=true)
    println()

    # check conservation
    tracer_total = PB.get_data(paleorun.output, "ocean.Tracer_total")
    @test isapprox(tracer_total[1], Tracer_total, rtol=1e-16)
    @test isapprox(tracer_total[end], tracer_total[1], rtol=1e-4)

    if do_plots
        pager = PALEOmodel.PlotPager((2, 2))
    else
        pager = PALEOmodel.NoPlotPager()
    end
   
    for t in [0.0, 0.5]
        pager(
            heatmap(paleorun.output, "oceansurface.surface_downwelling_photosynthetic_radiative_flux", (tmodel=t,), swap_xy=true),
            heatmap(paleorun.output, "oceansurface.open_area_fraction", (tmodel=t,), swap_xy=true),
            heatmap(paleorun.output, "ocean.temp", (tmodel=t, k=1), swap_xy=true),
        )     
    end

    for t in [0.0, 1.0]
        pager(
            heatmap(paleorun.output, "ocean.Tracer_conc", (tmodel=t, k=1), swap_xy=true),
        )
    end

    pager(
        plot(title="Tracer total", paleorun.output, ["ocean.Tracer_total"], ylabel="total (mol)",)
    )
  
end



end # testset
