using Test
using Logging
using DiffEqBase
using Sundials
import DataFrames

import PALEOboxes as PB

import PALEOocean
import PALEOmodel


@testset "romglb examples" begin

skipped_testsets = [
    # "O2_only", 
    # "P_O2",
    # "P_O2_S_Carb_open",
    # "P_O2_S_Carb_open_implicit", # requires PALEOaqchem > v0.3.9
]

matdir = joinpath(@__DIR__, "romaniello2010_transport") # assume zip file has been downloaded and unpacked in this subfolder

include("../atmreservoirreaction.jl")
include("SedimentationRate_dev.jl")

include("config_ocean_romglb_expts.jl")

!("O2_only" in skipped_testsets) && @testset "O2_only" begin

    model = PB.create_model_from_config(
        joinpath(@__DIR__, "PALEO_examples_romglb_cfg.yaml"), "romglb_abiotic_O2";
        modelpars=Dict(
            "matdir"=>matdir,
        )
    )
    tspan=(0,1e5)

    initial_state, modeldata = PALEOmodel.initialize!(model)
    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    PALEOmodel.ODE.integrateForwardDiff(
        paleorun, initial_state, modeldata, tspan, solvekwargs=(reltol=1e-9,)
    )

    println("conservation checks:")
    conschecks = [
        ("global", "total_O2",             1e-6 ),
    ]    
    for (domname, varname, rtol) in conschecks
        startval, endval = PB.get_data(paleorun.output, domname*"."*varname)[[1, end]]
        println("  check $domname.$varname $startval $endval $rtol")
        @test isapprox(startval, endval, rtol=rtol)
    end

    println("check values at end of run:")
    checkvals = [   
        ("atm", "pO2PAL",  0.995817,        1e-4),
        ("atm",  "O2_sms", 1000.0,          2.0),
    ]    
    for (domname, varname, checkval, rtol) in checkvals
        outputval = PB.get_data(paleorun.output, domname*"."*varname)[end]
        println("  check $domname.$varname $outputval $checkval $rtol")
        @test isapprox(outputval, checkval, rtol=rtol)
    end

end

!("P_O2" in skipped_testsets) && @testset "P_O2" begin

    model = PB.create_model_from_config(
        joinpath(@__DIR__, "PALEO_examples_romglb_cfg.yaml"), "romglb_P_O2";
        modelpars=Dict(
            "matdir"=>matdir,
        )
    )
    tspan=(0,1e5)

    initial_state, modeldata = PALEOmodel.initialize!(model)
    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    PALEOmodel.ODE.integrateForwardDiff(
        paleorun, initial_state, modeldata, tspan, solvekwargs=(reltol=1e-9,)
    )

    println("conservation checks:")
    conschecks = [
        ("global", "total_O2",             1e-6 ),
        ("ocean", "P_total",             1e-6 ),
    ]    
    for (domname, varname, rtol) in conschecks
        startval, endval = PB.get_data(paleorun.output, domname*"."*varname)[[1, end]]
        println("  check $domname.$varname $startval $endval $rtol")
        @test isapprox(startval, endval, rtol=rtol)
    end

    println("check values at end of run:")
    checkvals = [
        ("atm", "pO2PAL",  1.00155,                   1e-4),
        ("atm",  "O2_sms", 1000.0,                 2.0),
    ]    
    for (domname, varname, checkval, rtol) in checkvals
        outputval = PB.get_data(paleorun.output, domname*"."*varname)[end]
        println("  check $domname.$varname $outputval $checkval $rtol")
        @test isapprox(outputval, checkval, rtol=rtol)
    end

end

!("P_O2_S_Carb_open" in skipped_testsets) && @testset "P_O2_S_Carb_open" begin

    model = PB.create_model_from_config(
        joinpath(@__DIR__, "PALEO_examples_romglb_cfg.yaml"), "romglb_P_O2_S_Carb_open";
        modelpars=Dict(
            "matdir"=>joinpath(@__DIR__, "romaniello2010_transport"), # assume zip file has been downloaded and unpacked in this subfolder
            # Option 1: add pHfree, TAlk to state variables and apply constraint on TAlk      
            "TAlkStateExplicit"=>true,
            # Option 2: TAlk is an implicit variable, pHfree is a state variable
            # "TAlkStateExplicit"=>false,
        ),
    )

    # additional configuration to set shelf and deep carbonate burial
    config_ocean_romglb_expts(
        model, 
        [
            "set_carbonate_config_79box",  # configure shelf and deep carbonate burial for modern Earth
        ],
    )

    tspan=(0,1e5)

    initial_state, modeldata = PALEOmodel.initialize!(model)
    run = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    PALEOmodel.ODE.integrateDAEForwardDiff(
        run, initial_state, modeldata, tspan, 
        solvekwargs=(reltol=1e-6,)
        # solvekwargs=(reltol=1e-5,)
    )
    
    ncoutfile = tempname(cleanup=true)
    PALEOmodel.OutputWriters.save_netcdf(paleorun.output, ncoutfile; check_ext=false)
    ncoutput = PALEOmodel.OutputWriters.load_netcdf!(PALEOmodel.OutputWriters.OutputMemory(), ncoutfile; check_ext=false)

    for output in (paleorun.output, ncoutput)

        println("check values at end of run:")
        checkvals = [   
            ("fluxRtoOcean",    "flux_P",           4.17847e10,     1e-4),  # restoring flux
            ("atm",             "O2_restoring",    -2.6115e12,      1e-4),
            ("fluxOceanBurial", "flux_total_P",     4.17847e10,     1e-4),
            ("fluxOceanBurial", "flux_total_Corg",  2.61155e12,     1e-4),
            ("fluxOceanBurial", "flux_total_Ccarb", 2.3240e13,      1e-4),
        ]    
        for (domname, varname, checkval, rtol) in checkvals
            outputval = PB.get_data(output, domname*"."*varname)[end]
            println("  check $domname.$varname $outputval $checkval $rtol")
            @test isapprox(outputval, checkval, rtol=rtol)
        end
    end

end

!("P_O2_S_Carb_open_implicit" in skipped_testsets) && @testset "P_O2_S_Carb_open_implicit" begin

    model = PB.create_model_from_config(
        joinpath(@__DIR__, "PALEO_examples_romglb_cfg.yaml"), "romglb_P_O2_S_Carb_open";
        modelpars=Dict(
            "matdir"=>joinpath(@__DIR__, "romaniello2010_transport"), # assume zip file has been downloaded and unpacked in this subfolder
            # Option 1: add pHfree, TAlk to state variables and apply constraint on TAlk      
            # "TAlkStateExplicit"=>true,
            # Option 2: TAlk is an implicit variable, pHfree is a state variable
            "TAlkStateExplicit"=>false,
        ),
    )

    # additional configuration to set shelf and deep carbonate burial
    config_ocean_romglb_expts(
        model, 
        [
            "set_carbonate_config_79box",  # configure shelf and deep carbonate burial for modern Earth
        ],
    )

    tspan=(0,1e5)

    initial_state, modeldata = PALEOmodel.initialize!(model)
    run = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    PALEOmodel.ODE.integrateDAEForwardDiff(
        run, initial_state, modeldata, tspan, solvekwargs=(reltol=1e-5,)
    )

    println("check values at end of run:")
    checkvals = [ 
        ("fluxRtoOcean",    "flux_P",           4.17847e10,     1e-4),  # restoring flux
        ("atm",             "O2_restoring",    -2.6115e12,      1e-4),
        ("fluxOceanBurial", "flux_total_P",     4.17847e10,     1e-4),
        ("fluxOceanBurial", "flux_total_Corg",  2.61155e12,     1e-4),
        ("fluxOceanBurial", "flux_total_Ccarb", 2.3240e13,      1e-4),
    ]    
    for (domname, varname, checkval, rtol) in checkvals
        outputval = PB.get_data(run.output, domname*"."*varname)[end]
        println("  check $domname.$varname $outputval $checkval $rtol")
        @test isapprox(outputval, checkval, rtol=rtol)
    end

end


end
