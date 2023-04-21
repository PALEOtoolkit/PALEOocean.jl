
using Test
using Logging
using DiffEqBase
using Sundials
using Plots

import PALEOboxes as PB
import PALEOocean
import PALEOmodel

@testset "Ocean3box" begin

include("../examples/atmreservoirreaction.jl") # temporary solution to make ReactionReservoirAtm available

skipped_testsets = [
    # "Atm Ocean 3 box CO2", 
    # "Ocean 3 box remin O2",
    # "Ocean 3 box remin Ponly",
    # "Ocean 3 box bioprod P restore",
    # "Ocean 3 box airsea O2",
]

configfile = joinpath(@__DIR__, "configocean3box.yaml")

!("Atm Ocean 3 box CO2" in skipped_testsets) && @testset "Atm Ocean 3 box CO2" begin

    model = PB.create_model_from_config(configfile, "test_airsea_CO2")

    initial_state, modeldata = PALEOmodel.initialize!(model)

    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())
    # DAE with ForwardDiff sparse Jacobian

    PALEOmodel.ODE.integrateDAEForwardDiff(paleorun, initial_state, modeldata, (0, 1e5), alg=IDA(linear_solver=:KLU))
    

    for domainname in ("global", "atm", "ocean")
        println(domainname, " variables after integrate to steady-state:")
        domain = PB.get_domain(model, domainname)
        vars = PB.get_variables(domain)
        for var in vars
            println("\t", PB.fullname(var), " = ", PB.get_data(var, modeldata))
        end
    end

    # plot concentration
    display(plot(title="Total C", paleorun.output, ["global.total_C"], ylabel="C (mol)",))
    display(plot(title="Total C moldelta", paleorun.output, ["global.total_C.v_moldelta"], ylabel="C_moldelta (mol * delta)",))
    plot(title="Total C components", paleorun.output, ["global.total_C", "atm.CO2"], ylabel="C (mol)",)
    display(plot!(paleorun.output, "ocean.DIC", (cell=[:s, :h, :d], )))
    plot(title="Air-sea CO2", paleorun.output, "fluxAtmtoOceansurface.flux_CO2", (cell=[:s, :h], ), ylabel="mol yr-1", xlim=(0, 1e3))
    display(plot!(paleorun.output, "fluxAtmtoOceansurface.flux_total_CO2"))
    display(plot(title="pH (total)", paleorun.output, "ocean.pHtot", (cell=[:s, :h, :d], ), ylabel="pH (total)", xlim=(0, 1e3)))
    plot(title="d13C", paleorun.output, "atm.CO2_delta",ylabel="per mil", xlim=(0, 1e3))
    display(plot(paleorun.output, "ocean.DIC_delta", (cell=[:s, :h, :d], )))

    # test conservation
    total_C = PB.get_property(
        PB.get_data(paleorun.output, "global.total_C"),
        propertyname=:v
    )
    @test abs(total_C[end] - total_C[1])/total_C[1] < 1e-10
    total_C_moldelta = PB.get_property(
        PB.get_data(paleorun.output, "global.total_C"),
        propertyname=:v_moldelta
    )
    @test abs(total_C_moldelta[end] - total_C_moldelta[1])/total_C_moldelta[1] < 1e-10

    # test pCO2
    sigdigits = 10
    @test round.(PB.get_data(paleorun.output,  "atm.pCO2atm")[end], sigdigits=sigdigits) == 
        round(0.000651920453771681, sigdigits=sigdigits)
   
end


!("Ocean 3 box remin O2" in skipped_testsets) && @testset "Ocean 3 box remin O2" begin

    model = PB.create_model_from_config(configfile, "test_reminO2")

    ocean_domain = PB.get_domain(model, "ocean")
    @test PB.get_length(ocean_domain) == 3    

    initial_state, modeldata = PALEOmodel.initialize!(model)

    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    ocean_modelcreated_vars_dict = Dict([(var.name, var) for var in PB.get_variables(ocean_domain, hostdep=false)])
    
    
    # Check model derivative
    
    PB.do_deriv(modeldata.dispatchlists_all)

    println("state, sms variables after check model derivative:")
    for (state_var, sms_var) in PB.IteratorUtils.zipstrict(
                PB.get_vars(modeldata.solver_view_all.stateexplicit),
                PB.get_vars(modeldata.solver_view_all.stateexplicit_deriv)
            )
        println(PB.fullname(state_var), " ", PB.get_data(state_var, modeldata))
        println(PB.fullname(sms_var), " ", PB.get_data(sms_var, modeldata))
    end

    println("ocean  model created variables after check model derivative:")
    for (name, var ) in ocean_modelcreated_vars_dict
        println("\t", PB.fullname(var), " = ", PB.get_data(var, modeldata))
    end

    # check conservation
    ocean_P_sms_data = PB.get_data(PB.get_variable(ocean_domain,"P_sms"), modeldata)
    sum_P_sms = sum(ocean_P_sms_data)
    println("check conservation: sum(P_sms)=",sum_P_sms)
    @test abs(sum_P_sms) < 1e-16

    ocean_O2_sms_data = PB.get_data(PB.get_variable(ocean_domain,"O2_sms"), modeldata)
    sum_O2_sms = sum(ocean_O2_sms_data)
    println("check conservation: sum(O2_sms)=",sum_O2_sms)
    @test abs(sum_O2_sms) < 1e-16

    # integrate to approx steady state
    PALEOmodel.ODE.integrate(paleorun, initial_state, modeldata, (0, 1e5))  # first run includes JIT time
    
    # show(paleorun.output.domains["ocean"].data, allcols=true)
    # println()

    # check conservation
    P_total = PB.get_data(paleorun.output, "ocean.P_total")
    @test abs(P_total[end] - P_total[1])/P_total[1] < 1e-10

    O2_total = PB.get_data(paleorun.output, "ocean.O2_total")
    @test abs(O2_total[end] - O2_total[1])/O2_total[1] < 1e-10

    DIC_total = PB.get_property(
        PB.get_data(paleorun.output, "ocean.DIC_total"),
        propertyname=:v
    )
    @test abs(DIC_total[end] - DIC_total[1])/DIC_total[1] < 1e-10
    DIC_total_moldelta = PB.get_property(
        PB.get_data(paleorun.output, "ocean.DIC_total"),
        propertyname=:v_moldelta
    )
    @test abs(DIC_total_moldelta[end] - DIC_total_moldelta[1])/DIC_total_moldelta[1] < 1e-10

    TAlk_total = PB.get_data(paleorun.output, "ocean.TAlk_total")
    @test abs(TAlk_total[end] - TAlk_total[1])/TAlk_total[1] < 1e-10
end



!("Ocean 3 box remin Ponly" in skipped_testsets) && @testset "Ocean 3 box remin Ponly" begin

    model = PB.create_model_from_config(configfile, "test_reminPonly")


    ocean_domain = PB.get_domain(model, "ocean")
    @test PB.get_length(ocean_domain) == 3

    initial_state, modeldata = PALEOmodel.initialize!(model)

    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    ocean_modelcreated_vars_dict = Dict([(var.name, var) for var in PB.get_variables(ocean_domain, hostdep=false)])
    
    
    # Check model derivative
    
    PB.do_deriv(modeldata.dispatchlists_all)

    println("state, sms variables after check model derivative:")
    for (state_var, sms_var) in PB.IteratorUtils.zipstrict(
                PB.get_vars(modeldata.solver_view_all.stateexplicit),
                PB.get_vars(modeldata.solver_view_all.stateexplicit_deriv)
            )
        println(PB.fullname(state_var), " ", PB.get_data(state_var, modeldata))
        println(PB.fullname(sms_var), " ", PB.get_data(sms_var, modeldata))
    end

    println("ocean  model created variables after check model derivative:")
    for (name, var ) in ocean_modelcreated_vars_dict
        println("\t", PB.fullname(var), " = ", PB.get_data(var, modeldata))
    end

    # check conservation
    ocean_P_sms_data = PB.get_data(PB.get_variable(ocean_domain,"P_sms"), modeldata)
    sum_P_sms = sum(ocean_P_sms_data)
    println("check conservation: sum(P_sms)=",sum_P_sms)
    @test abs(sum_P_sms) < 1e-16

    # integrate to approx steady state
    PALEOmodel.ODE.integrate(paleorun, initial_state, modeldata, (0, 1e5))  # first run includes JIT time
    
    # show(paleorun.output["ocean"], allcols=true)
    println()

    # check conservation
    P_total = PB.get_data(paleorun.output, "ocean.P_total")
    @test abs(P_total[end] - P_total[1])/P_total[1] < 1e-10
end



!("Ocean 3 box bioprod P restore" in skipped_testsets) && @testset "Ocean 3 box bioprod P restore" begin

    model = PB.create_model_from_config(configfile, "test_bioprodPrest")

    ocean_domain = PB.get_domain(model, "ocean")
    @test PB.get_length(ocean_domain) == 3

    initial_state, modeldata = PALEOmodel.initialize!(model)

    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    ocean_modelcreated_vars_dict = Dict([(var.name, var) for var in PB.get_variables(ocean_domain, hostdep=false)])
    
    
    # Check model derivative
    
    PB.do_deriv(modeldata.dispatchlists_all)

    println("state, sms variables after check model derivative:")
    for (state_var, sms_var) in PB.IteratorUtils.zipstrict(
                PB.get_vars(modeldata.solver_view_all.stateexplicit),
                PB.get_vars(modeldata.solver_view_all.stateexplicit_deriv)
            )
        println(PB.fullname(state_var), " ", PB.get_data(state_var, modeldata))
        println(PB.fullname(sms_var), " ", PB.get_data(sms_var, modeldata))
    end

    println("ocean  model created variables after check model derivative:")
    for (name, var ) in ocean_modelcreated_vars_dict
        println("\t", PB.fullname(var), " = ", PB.get_data(var, modeldata))
    end

    # check conservation
    ocean_P_sms_data = PB.get_data(PB.get_variable(ocean_domain,"P_sms"), modeldata)
    sum_P_sms = sum(ocean_P_sms_data)
    println("check conservation: sum(P_sms)=",sum_P_sms)
    @test abs(sum_P_sms) < 1e-16

    # integrate to approx steady state
    PALEOmodel.ODE.integrate(paleorun, initial_state, modeldata, (0, 1e5))  # first run includes JIT time
    
    # show(paleorun.output["ocean"], allcols=true)
    # println()

    # check conservation
    P_total = PB.get_data(paleorun.output, "ocean.P_total")
    @test abs(P_total[end] - P_total[1])/P_total[1] < 1e-10
end



!("Ocean 3 box airsea O2" in skipped_testsets) && @testset "Ocean 3 box airsea O2" begin

    model = PB.create_model_from_config(configfile, "test_airsea_O2")

    # Test OceanBase domain configuration
    @test PB.get_num_domains(model) == 6

    global_domain = PB.get_domain(model, "global")
    @test PB.get_length(global_domain) == 1

    ocean_domain = PB.get_domain(model, "ocean")
    @test PB.get_length(ocean_domain) == 3

    @test ocean_domain.grid.subdomains["oceansurface"].indices == [1, 2]
    @test ocean_domain.grid.subdomains["oceanfloor"].indices == [1, 2, 3]

    oceansurface_domain = PB.get_domain(model, "oceansurface")
    @test PB.get_length(oceansurface_domain) == 2

    oceanfloor_domain = PB.get_domain(model, "oceanfloor")
    @test PB.get_length(oceanfloor_domain) == 3
    
    # test OceanBase variables
    
    initial_state, modeldata = PALEOmodel.initialize!(model)

    paleorun = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    ocean_modelcreated_vars_dict = Dict([(var.name, var) for var in PB.get_variables(ocean_domain, hostdep=false)])
    
    println("ocean  model created variables after initialize!:")
    for (name, var ) in ocean_modelcreated_vars_dict
        println("\t", PB.fullname(var), " = ", PB.get_data(var, modeldata))
    end

    # test a constant variable
    @test PB.get_data(ocean_modelcreated_vars_dict["zupper"], modeldata) == [0.0, 0.0, -100.0]

   
    # bodge a test for ocean tracer with single non-zero cell
    ocean_Tracer_data = PB.get_data(PB.get_variable(ocean_domain,"Tracer"), modeldata)
    ocean_Tracer_data .= [1.0, 0.0, 0.0]
    # update initial_state with our bodged values
    initial_state = PALEOmodel.get_statevar(modeldata.solver_view_all)

    # Check model derivative
    
    PB.do_deriv(modeldata.dispatchlists_all)

    println("state, sms variables after check model derivative:")
    for (state_var, sms_var) in PB.IteratorUtils.zipstrict(
                PB.get_vars(modeldata.solver_view_all.stateexplicit),
                PB.get_vars(modeldata.solver_view_all.stateexplicit_deriv)
            )
        println(PB.fullname(state_var), " ", PB.get_data(state_var, modeldata))
        println(PB.fullname(sms_var), " ", PB.get_data(sms_var, modeldata))
    end

    # check conservation
    ocean_Tracer_sms_data = PB.get_data(PB.get_variable(ocean_domain,"Tracer_sms"), modeldata)
    sum_Tracer_sms = sum(ocean_Tracer_sms_data)
    println("check conservation: sum(Tracer_sms)=",sum_Tracer_sms)
    @test abs(sum_Tracer_sms) < 1e-16

    # integrate to approx steady state
    PALEOmodel.ODE.integrate(paleorun, initial_state, modeldata, (0, 1e5))  # first run includes JIT time
    #show(paleorun.output["fluxAtmtoOceansurface"], allcols=true)
    #println()
    #show(paleorun.output["atm"], allcols=true)
    #println()
    #show(paleorun.output["ocean"], allcols=true)
    #println()

    # check conservation
    tracer_total = PB.get_data(paleorun.output, "ocean.Tracer_total")
    @test abs(tracer_total[1] - 1.0) < 1e-16
    @test abs(tracer_total[end] - tracer_total[1]) < 1e-4

    total_O2 = PB.get_data(paleorun.output, "global.total_O2")
    @test abs(total_O2[end] - total_O2[1]) < 1e-7*total_O2[1]

    # plot concentration
    display(plot(title="Tracer total", paleorun.output, ["ocean.Tracer_total"], ylabel="total (mol)",))
    display(plot(title="Tracer concentration", paleorun.output, "ocean.Tracer_conc", (cell=[:s, :h, :d], ), ylabel="conc (mol m-3)", xlim=(0, 1e3)))
    display(plot(title="Total O2", paleorun.output, ["global.total_O2"], ylabel="O2 (mol)",))
    plot(title="Total O2 components", paleorun.output, ["global.total_O2", "atm.O"], ylabel="O2 (mol)",)
    display(plot!(paleorun.output, "ocean.O2", (cell=[:s, :h, :d], )))
    plot(title="Air-sea O2", paleorun.output, "fluxAtmtoOceansurface.flux_O2", (cell=[:s, :h], ), ylabel="mol yr-1", xlim=(0, 1e3))
    display(plot!(paleorun.output, "fluxAtmtoOceansurface.flux_total_O2"))

    model = nothing
end

end

