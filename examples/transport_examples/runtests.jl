using Logging
using Test

import PALEOboxes as PB
import PALEOmodel
import PALEOocean

include("TransportExamples.jl")

@testset "transport_examples" begin

skipped_testsets = [
    # "transport_advect",
    # "transport_diffuse",
]

!("transport_advect" in skipped_testsets) && @testset "transport_advect" begin

    model = PB.create_model_from_config(
        joinpath(@__DIR__, "TransportExamples_cfg.yaml"), "example_advect"
    )

    initial_state, modeldata = PALEOmodel.initialize!(model; check_units_opt=:error)

    # bodge an updated initial_state for testing
    ocean_T = PB.get_data(PB.get_variable(model, "ocean.T"), modeldata) # model data array for tracer T
    ocean_T[1] = 2e14*200.0*1.0  # ~ 1 mol m-3 in first cell (top of first column)
    initial_state = PALEOmodel.get_statevar(modeldata.solver_view_all)

    tspan=(0.0, 5e3)

    run = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    PALEOmodel.ODE.integrateForwardDiff(
        run, initial_state, modeldata, tspan,
        solvekwargs=(reltol=1e-5,),
    )

    println("conservation checks:")
    conschecks = [
        ("ocean",  "T_total", 1e-14),
    ]
    for (domname, varname, rtol) in conschecks
        startval, endval = PB.get_data(run.output, domname*"."*varname)[[1, end]]
        println("  check $domname.$varname $startval $endval $rtol")
        @test isapprox(startval, endval, rtol=rtol)
    end

end

!("transport_diffuse" in skipped_testsets) && @testset "transport_diffuse" begin

    model = PB.create_model_from_config(
        joinpath(@__DIR__, "TransportExamples_cfg.yaml"), "example_diffuse"
    )

    initial_state, modeldata = PALEOmodel.initialize!(model; check_units_opt=:error)

    # bodge an updated initial_state for testing
    ocean_T = PB.get_data(PB.get_variable(model, "ocean.T"), modeldata) # model data array for tracer T
    ocean_T[50] = 2e14*20.0*1.0  # ~ 1 mol m-3 in middle cell of first column
    ocean_T[100+5] = 2e13*200.0*1.0  # ~ 1 mol m-3 in middle cell of second column
    initial_state = PALEOmodel.get_statevar(modeldata.solver_view_all)

    tspan=(0.0, 1e6)

    run = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    PALEOmodel.ODE.integrateForwardDiff(
        run, initial_state, modeldata, tspan,
        solvekwargs=(reltol=1e-5,),
    )

    println("conservation checks:")
    conschecks = [
        ("ocean",  "T_total", 1e-14),
    ]
    for (domname, varname, rtol) in conschecks
        startval, endval = PB.get_data(run.output, domname*"."*varname)[[1, end]]
        println("  check $domname.$varname $startval $endval $rtol")
        @test isapprox(startval, endval, rtol=rtol)
    end

    println("check values at end of run:")
    checkvals = [  
        ("ocean", "T_conc",  0.01*ones(110),    1e-5 ),  # check uniform concentration in 110 ocean cells
    ]    
    for (domname, varname, checkval, rtol) in checkvals
        outputval = PB.get_data(run.output, domname*"."*varname)[end]
        println("  check $domname.$varname $outputval $checkval $rtol")
        @test isapprox(outputval, checkval, rtol=rtol)
    end

end



end
