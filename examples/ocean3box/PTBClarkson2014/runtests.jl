using Test
using Logging
using DiffEqBase
using Sundials
import DataFrames

import PALEOboxes as PB

import PALEOreactions
import PALEOmodel


@testset "PTB examples" begin

skipped_testsets = [
    # "PTB_2015",
]

!("PTB_2015" in skipped_testsets) && @testset "PTB_2014" begin

    include("config_PTB3box_expts.jl")

    # Clarkson etal (2015) CO2LO scenario
    model = config_PTB3box_expts("Co2LOmHWC4pp", ["Sw_2Ts", "Pp_PEes", "Lk_2", "Cia_s2", "Cib_1"])
    tspan=(-260e6, -251.88e6) # stop just after second C pulse # (-260e6,-240e6) 

    initial_state, modeldata = PALEOmodel.initialize!(model)
    run = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

    sol = PALEOmodel.ODE.integrateDAEForwardDiff(
        run, initial_state, modeldata, tspan,
        alg=IDA(linear_solver=:KLU),
        solvekwargs=(
            abstol=1e-6*PALEOmodel.get_statevar_norm(modeldata.solver_view_all),
            reltol=1e-5,
            save_start=false,
            dtmax=0.5e5,
        )
    )

    println("conservation checks:")
    conschecks = [ 
        ("global", "total_S",       :v,      1e-6 ),
        # ("global", "total_S",     :v_moldelta,    1e-6), # starts at 0.0 so rtol isn't useful
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
        ("ocean",    "pHtot",        [7.4870, 7.4438, 7.1207],  1e-4 ),
        ("atm",     "pCO2atm",      5926.39e-6,                 2.5e-3),
        ("atm",     "CO2_delta",    -7.0371,                    1e-4),
        ("ocean",   "P_total",      6.8448e15,                  1e-3),
    ]    
    for (domname, varname, checkval, rtol) in checkvals
        outputval = PB.get_data(run.output, domname*"."*varname)[end]
        println("  check $domname.$varname $outputval $checkval $rtol")
        @test isapprox(outputval, checkval, rtol=rtol)
    end

end


end
