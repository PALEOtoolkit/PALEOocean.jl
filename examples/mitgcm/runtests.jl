using Test
using Logging
import DataFrames

import PALEOboxes as PB

import PALEOocean
import PALEOmodel


@testset "MITgcm examples" begin

skipped_testsets = [
    # "PO4MMbase",
]


!("PO4MMbase" in skipped_testsets) && @testset "PO4MMbase" begin

    include("config_mitgcm_expts.jl")

    model = PB.create_model_from_config(
        joinpath(@__DIR__, "MITgcm_2deg8_PO4MMbase.yaml"), "PO4MMbase"
    )
    toutputs = [0.0, 0.25, 0.5, 0.75, 1.0]

    tstep_explicit_s::Int = 86400 # s timestep to use for explicit transport
    tstep_implicit_s::Int = tstep_explicit_s
    tstep_explicit_yr = tstep_explicit_s/PB.Constants.k_secpyr # yr
    @info "using timesteps tstep_implicit_s $tstep_implicit_s tstep_explicit_s $tstep_explicit_s tstep_explicit_yr $tstep_explicit_yr"
    transportMITgcm = PB.get_reaction(model, "ocean", "transportMITgcm")
    PB.setvalue!(transportMITgcm.pars.Aimp_deltat, tstep_implicit_s)
    transportMITgcm = nothing # holds large transport matrix arrays !


    initial_state, modeldata = PALEOmodel.initialize!(model)
    
    paleorun = PALEOmodel.Run(model=model, output=PALEOmodel.OutputWriters.OutputMemory())

    @time PALEOmodel.ODEfixed.integrateEuler(paleorun, initial_state, modeldata, toutputs, tstep_explicit_yr)

    println("conservation checks:")
    conschecks = [
        ("global", "total_O2",          1e-6),
        ("global",  "total_P",          1e-6),
    ]    
    for (domname, varname, rtol) in conschecks
        startval, endval = PB.get_data(paleorun.output, domname*"."*varname)[[1, end]]
        println("  check $domname.$varname $startval $endval $rtol")
        @test isapprox(startval, endval, rtol=rtol)
    end

    println("check values at end of run:")
    checkvals = [
        ("ocean", "Prod_Corg_total",  3.6975e15,    1e-4),
    ]    
    for (domname, varname, checkval, rtol) in checkvals
        outputval = PB.get_data(paleorun.output, domname*"."*varname)[end]
        println("  check $domname.$varname $outputval $checkval $rtol")
        @test isapprox(outputval, checkval, rtol=rtol)
    end

end



end
