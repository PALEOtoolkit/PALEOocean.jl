using Test
using Logging
using DiffEqBase
using Sundials
import DataFrames
import SparseArrays
import SparseDiffTools

import PALEOboxes as PB

import PALEOocean
import PALEOmodel
import PALEOcopse


@testset "romglb examples" begin

skipped_testsets = [
    # "airsea_O2",
    # "O2_only", 
    # "P_O2",
    # "P_O2_S_Carb_open",
    # "P_O2_S_Carb_open_implicit", # requires PALEOaqchem > v0.3.9
]

matdir = joinpath(@__DIR__, "romaniello2010_transport") # assume zip file has been downloaded and unpacked in this subfolder

if !isdir(matdir)
    include("download_romaniello2010_files.jl")
    download_romaniello2010_files()
end

include("../atmreservoirreaction.jl")
include("SedimentationRate_dev.jl")

include("config_ocean_romglb_expts.jl")

!("airsea_O2" in skipped_testsets) && @testset "airsea_O2" begin

    model = PB.create_model_from_config(
        joinpath(@__DIR__, "PALEO_examples_romglb_cfg.yaml"), "romglb_abiotic_O2";
        modelpars=Dict(
            "matdir"=>matdir,
        )
    )

    # Test OceanBase domain configuration
    @test PB.get_num_domains(model) == 6

    global_domain = PB.get_domain(model, "global")
    @test PB.get_length(global_domain) == 1

    ocean_domain = PB.get_domain(model, "ocean")
    @test PB.get_length(ocean_domain) == 79

    
    # test OceanBase variables

    initial_state, modeldata = PALEOmodel.initialize!(model)

    ocean_modelcreated_vars_dict = Dict([(var.name, var) for var in PB.get_variables(ocean_domain, hostdep=false)])
    
    println("ocean  model created variables after initialize!:")
    for (name, var ) in ocean_modelcreated_vars_dict
        println("\t", PB.fullname(var), " = ", PB.get_data(var, modeldata))
    end

    # bodge a test for ocean tracer with single non-zero cell
    ocean_T_data = PB.get_data(PB.get_variable(ocean_domain,"T"), modeldata)
    ocean_T_data .= 0.0
    ocean_T_data[1] = 1.0
    # update initial_state with our bodged values
    initial_state = PALEOmodel.get_statevar(modeldata.solver_view_all)
     
    # Check model derivative
    
    PB.do_deriv(modeldata.dispatchlists_all)

    println("state, sms variables after check model derivative:")
    for (state_var, sms_var) in PB.IteratorUtils.zipstrict(PB.get_vars(modeldata.solver_view_all.stateexplicit), 
                                    PB.get_vars(modeldata.solver_view_all.stateexplicit_deriv))
        println(PB.fullname(state_var), " ", PB.get_data(state_var, modeldata))
        println(PB.fullname(sms_var), " ", PB.get_data(sms_var, modeldata))
    end

    # check conservation
    ocean_T_sms_data = PB.get_data(PB.get_variable(ocean_domain,"T_sms"), modeldata)
    sum_T_sms = sum(ocean_T_sms_data)
    println("check conservation: sum(Tracer_sms)=",sum_T_sms)
    @test abs(sum_T_sms) < 1e-15

    # check jacobian sparsity calculation
    jac_prototype = PALEOmodel.JacobianAD.calcJacobianSparsitySparsityTracing!(model, modeldata, initial_state, 0.0)
    @test SparseArrays.nnz(jac_prototype) == 819
    jac_proto_fill = PALEOmodel.SparseUtils.fill_sparse_jac(jac_prototype)
    colors = SparseDiffTools.matrix_colors(copy(jac_proto_fill))
    @test maximum(colors) == 23  

    # integrate to approx steady state
    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())
    PALEOmodel.ODE.integrateForwardDiff(paleorun, initial_state, modeldata, (0, 1e5), solvekwargs=(reltol=1e-5,))  # first run includes JIT time
    
    # check conservation
    T_total = PB.get_data(paleorun.output, "ocean.T_total")
    @test abs(T_total[1] - 1.0) < 1e-16
    @test abs(T_total[end] - T_total[1]) < 1e-4

    total_O2 = PB.get_data(paleorun.output, "global.total_O2")
    @test abs(total_O2[end] - total_O2[1]) < 1e-7*total_O2[1]

end

!("O2_only" in skipped_testsets) && @testset "O2_only" begin

    model = PB.create_model_from_config(
        joinpath(@__DIR__, "PALEO_examples_romglb_cfg.yaml"), "romglb_abiotic_O2";
        modelpars=Dict(
            "matdir"=>matdir,
        )
    )
    tspan=(0,1e5)

    initial_state, modeldata = PALEOmodel.initialize!(model; check_units_opt=:error)
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

    initial_state, modeldata = PALEOmodel.initialize!(model; check_units_opt=:error)
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
            "matdir"=>matdir, # assume zip file has been downloaded and unpacked in this subfolder
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

    initial_state, modeldata = PALEOmodel.initialize!(model; check_units_opt=:error)
    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    PALEOmodel.ODE.integrateDAEForwardDiff(
        paleorun, initial_state, modeldata, tspan, 
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
            "matdir"=>matdir, # assume zip file has been downloaded and unpacked in this subfolder
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

    initial_state, modeldata = PALEOmodel.initialize!(model; check_units_opt=:error)
    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    PALEOmodel.ODE.integrateDAEForwardDiff(
        paleorun, initial_state, modeldata, tspan, solvekwargs=(reltol=1e-5,)
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
        outputval = PB.get_data(paleorun.output, domname*"."*varname)[end]
        println("  check $domname.$varname $outputval $checkval $rtol")
        @test isapprox(outputval, checkval, rtol=rtol)
    end

end


end

nothing # so no output printed when run from REPL