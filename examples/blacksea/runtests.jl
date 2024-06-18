using Test
using Logging
using DiffEqBase
using Sundials
import DataFrames

import PALEOboxes as PB

import PALEOocean
import PALEOmodel


@testset "rom Black_Sea examples" begin

skipped_testsets = [
    # "O2_only", 
    # "P_O2",
    
]

matdir = joinpath(@__DIR__, "../romglb/romaniello2010_transport") # assume zip file has been downloaded and unpacked in this subfolder

include("config_ocean_blacksea_expts.jl")


!("O2_only" in skipped_testsets) && @testset "O2_only" begin

    model = PB.create_model_from_config(
        joinpath(@__DIR__, "PALEO_examples_blacksea_cfg.yaml"), "blacksea_abiotic_O2",
        modelpars=Dict(
            "matdir"=>matdir,
        )
    )

    tspan=(0, 1e5)

    initial_state, modeldata = PALEOmodel.initialize!(model)
    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    PALEOmodel.ODE.integrateForwardDiff(
        paleorun, initial_state, modeldata, tspan, solvekwargs=(reltol=1e-9,)
    )

    println("conservation checks:")
    conschecks = [
        ("ocean", "T_total",             1e-6 ),
    ]    
    for (domname, varname, rtol) in conschecks
        startval, endval = PB.get_data(paleorun.output, domname*"."*varname)[[1, end]]
        println("  check $domname.$varname $startval $endval $rtol")
        @test isapprox(startval, endval, rtol=rtol)
    end

    println("check values at end of run:")
    checkvals = [   
    ]    
    for (domname, varname, checkval, rtol) in checkvals
        outputval = PB.get_data(paleorun.output, domname*"."*varname)[end]
        println("  check $domname.$varname $outputval $checkval $rtol")
        @test isapprox(outputval, checkval, rtol=rtol)
    end

end

!("P_O2" in skipped_testsets) && @testset "P_O2" begin

    model = PB.create_model_from_config(
        joinpath(@__DIR__, "PALEO_examples_blacksea_cfg.yaml"), "blacksea_P_O2",
        modelpars=Dict(
            "matdir"=>matdir,
        )
    )

    tspan=(0, 1e5)

    initial_state, modeldata = PALEOmodel.initialize!(model)
    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    PALEOmodel.ODE.integrateForwardDiff(
        paleorun, initial_state, modeldata, tspan, solvekwargs=(reltol=1e-9,)
    )

    ncoutfile = tempname(cleanup=true)
    PALEOmodel.OutputWriters.save_netcdf(paleorun.output, ncoutfile; check_ext=false)
    ncoutput = PALEOmodel.OutputWriters.load_netcdf!(PALEOmodel.OutputWriters.OutputMemory(), ncoutfile; check_ext=false)

    for output in (paleorun.output, ncoutput)
        println("conservation checks:")
        conschecks = [
            ("ocean", "T_total",             1e-6 ),
            ("ocean", "P_total",             1e-6 ),
        ]    
        for (domname, varname, rtol) in conschecks
            startval, endval = PB.get_data(output, domname*"."*varname)[[1, end]]
            println("  check $domname.$varname $startval $endval $rtol")
            @test isapprox(startval, endval, rtol=rtol)
        end

        println("check values at end of run:")
        checkvals = [
        ]    
        for (domname, varname, checkval, rtol) in checkvals
            outputval = PB.get_data(output, domname*"."*varname)[end]
            println("  check $domname.$varname $outputval $checkval $rtol")
            @test isapprox(outputval, checkval, rtol=rtol)
        end
    end

end


end
