using Test
using Logging
using DiffEqBase
using Sundials
import DataFrames

import PALEOboxes as PB

import PALEOocean
import PALEOcopse
import PALEOmodel


@testset "ocean3box examples" begin

skipped_testsets = [
    # "oaonly_abiotic",
    # "oaonly",
    # "oaopencarb",
]

!("oaonly_abiotic" in skipped_testsets) && @testset "oaonly_abiotic" begin

    include("config_ocean3box_expts.jl")

    model = config_ocean3box_expts("oaonly_abiotic", ["baseline"]); tspan=(0,1e4)

    initial_state, modeldata = PALEOmodel.initialize!(model)
    run = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

    sol = PALEOmodel.ODE.integrateDAEForwardDiff(
        run, initial_state, modeldata, tspan,
        alg=IDA(linear_solver=:KLU),
        solvekwargs=(
            abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all),
            save_start=false
        )
    )

    println("conservation checks:")
    conschecks = [
        ("global", "total_C", :v,           1e-6 ),
        ("global", "total_C", :v_moldelta,  1e-6),
        ("ocean",  "TAlk_total", nothing,   1e-6),
    ]
    for (domname, varname, propertyname, rtol) in conschecks
        startval, endval = PB.get_property(
            PB.get_data(run.output, domname*"."*varname),
            propertyname=propertyname
        )[[1, end]]
        println("  check $domname.$varname $startval $endval $rtol")
        @test isapprox(startval, endval, rtol=rtol)
    end

    println("check values at end of run:")
    checkvals = [  
        ("ocean", "pHfree",  [8.236714539086925, 8.22352510474905, 8.149122544874377],    1e-4 ),
        ("atm", "pCO2atm",  301.244e-6,                                                  1e-4),
        ("atm",  "CO2_delta", -10.161,                                                  1e-4),
    ]    
    for (domname, varname, checkval, rtol) in checkvals
        outputval = PB.get_data(run.output, domname*"."*varname)[end]
        println("  check $domname.$varname $outputval $checkval $rtol")
        @test isapprox(outputval, checkval, rtol=rtol)
    end

end

!("oaonly" in skipped_testsets) && @testset "oaonly" begin

    include("config_ocean3box_expts.jl")

    model = config_ocean3box_expts("oaonly", ["baseline"]); tspan=(0,1e4)

    initial_state, modeldata = PALEOmodel.initialize!(model)
    run = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

    sol = PALEOmodel.ODE.integrateDAEForwardDiff(
        run, initial_state, modeldata, tspan,
        alg=IDA(linear_solver=:KLU),
        solvekwargs=(
            abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all),
            save_start=false
        )
    )

    println("conservation checks:")
    conschecks = [  
        ("global", "total_C",   :v,             1e-6 ),
        ("global", "total_C",   :v_moldelta,    1e-6),
        ("global", "total_O2",  nothing,        1e-6 ),
        ("global", "total_S",   :v,             1e-6 ),
        ("global", "total_S",   :v_moldelta,    1.0 ), # starts at 0.0 so no reltol test
        ("ocean",  "TAlk_total",nothing,        1e-6)
    ]    
    for (domname, varname, propertyname, rtol) in conschecks
        startval, endval = PB.get_property(
            PB.get_data(run.output, domname*"."*varname),
            propertyname=propertyname
        )[[1, end]]
        println("  check $domname.$varname $startval $endval $rtol")
        @test isapprox(startval, endval, rtol=rtol)
    end

    println("check values at end of run:")
    checkvals = [
        ("ocean", "pHfree",  [8.38775021448627, 8.367024330312566, 8.107940989697116],   1e-4 ),
        ("atm", "pCO2atm",  191.423e-6,                                                 1e-4),
        ("atm",  "CO2_delta", -8.1965,                                                  1e-4),
        ("atm", "pO2atm",   0.2102,                                                     1e-4),
    ]    
    for (domname, varname, checkval, rtol) in checkvals
        outputval = PB.get_data(run.output, domname*"."*varname)[end]
        println("  check $domname.$varname $outputval $checkval $rtol")
        @test isapprox(outputval, checkval, rtol=rtol)
    end

end

!("oaopencarb" in skipped_testsets) && @testset "oaopencarb" begin

    include("config_ocean3box_expts.jl")

    model = config_ocean3box_expts("oaopencarb", ["killbio", "lowO2"]); tspan=(-10e6,10e6) 

    initial_state, modeldata = PALEOmodel.initialize!(model)
    run = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

    sol = PALEOmodel.ODE.integrateDAEForwardDiff(
        run, initial_state, modeldata, tspan,
        alg=IDA(linear_solver=:KLU),
        solvekwargs=(
            abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all),
            save_start=false
        )
    )

    println("conservation checks:")
    conschecks = [  
        ("global", "total_S",   :v,             1e-6 ),
        ("global", "total_S",   :v_moldelta,    1),  # run starts with 0.0
        ("ocean",  "P_total",   nothing,        1e-6)
    ]    
    for (domname, varname, propertyname, rtol) in conschecks
        startval, endval = PB.get_property(
            PB.get_data(run.output, domname*"."*varname),
            propertyname=propertyname
        )[[1, end]]
        println("  check $domname.$varname $startval $endval $rtol")
        @test isapprox(startval, endval, rtol=rtol)
    end

    println("check values at end of run:")
    checkvals = [ 
        ("ocean", "pHtot", [8.32404, 8.40819, 8.3423],    1e-4 ),
        ("atm", "pCO2atm",  0.000209573,            1e-4),
        ("atm", "CO2_delta", -8.51631 ,          1e-4),             
    ]    
    for (domname, varname, checkval, rtol) in checkvals
        outputval = PB.get_data(run.output, domname*"."*varname)[end]
        println("  check $domname.$varname $outputval $checkval $rtol")
        @test isapprox(outputval, checkval, rtol=rtol)
    end

end


end
