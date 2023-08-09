module Burial


import PALEOboxes as PB
using PALEOboxes.DocStrings

using SpecialFunctions
using Interpolations

# import Infiltrator

"""
    ReactionShelfCarb

Shallow-water carbonate burial controlled by carbonate saturation state (after [Caldeira1993](@cite))

Carbonate burial rate in cell `i` is `shelf_Ccarb[i]` = `carbsedshallow` * `shelfareanorm[i]` * '(OmegaAR[i]-1.0)^1.7`
where `carbsedshallow` (mol C yr-1) controls the global rate, and `shelfareanorm[i]` controls the spatial distribution 
among ocean shelf cells.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionShelfCarb{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("carbsedshallow", 1.4355e12, units="mol C yr-1",
            description="total carbonate deposition rate"),

        PB.ParDoubleVec("shelfareanorm", Float64[], units="",
            description="per box distribution of carbonate burial (length=Domain size, must sum to 1.0)"),

        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
    )

end

function PB.register_methods!(rj::ReactionShelfCarb)
    
    CIsotopeType = rj.pars.CIsotope[]
    PB.setfrozen!(rj.pars.CIsotope)
    @info "register_methods! $(PB.fullname(rj)) CIsotopeType=$CIsotopeType"

    fluxOceanfloorSolute = PB.Fluxes.FluxContrib(
        "fluxOceanfloor.soluteflux_", ["DIC::$CIsotopeType", "TAlk", "(Ca)"],
        alloptional=false
    )

    fluxOceanBurial = PB.Fluxes.FluxContrib(
        "fluxOceanBurial.flux_", ["Ccarb::$CIsotopeType"],
    )

    vars = [
        PB.VarDepScalar("(global.shelfarea_force)",   "",  
            "optional forcing multiplier for carbsedshallow (defaults to 1.0)"),       

        PB.VarDep("ocean.oceanfloor.OmegaAR", "", "aragonite saturation"),

        PB.VarProp("shelf_Ccarb", "mol yr-1", "shelf Ccarb burial",
            attributes=(:field_data=>CIsotopeType, :initialize_to_zero=>true, :calc_total=>true)),
    ]
    if CIsotopeType <: PB.AbstractIsotopeScalar
        push!(vars,
            PB.VarDep("ocean.oceanfloor.DIC_delta",          "per mil",  "d13C DIC"),
            PB.VarDepScalar("ocean.D_mccb_DIC",   "per mil",  
                "d13C marine calcite burial relative to ocean DIC"),
        )
    end

    PB.add_method_do!(
        rj,
        do_shelf_carb,
        (
            PB.VarList_namedtuple_fields(fluxOceanfloorSolute),
            PB.VarList_namedtuple_fields(fluxOceanBurial),
            PB.VarList_namedtuple(vars)
        ),
        p = CIsotopeType, 
    )

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end
    
function PB.check_configuration(rj::ReactionShelfCarb, model::PB.Model)
    configok = true

    PB.check_parameter_sum(rj.pars.shelfareanorm, PB.get_length(rj.domain)) || (configok = false)

    return configok
end

function do_shelf_carb(
    m::PB.ReactionMethod,
    pars,
    (fluxOceanfloorSolute, fluxOceanBurial, vars), 
    cellrange::PB.AbstractCellRange,
    delta
)
    CIsotopeType = m.p

    shelfarea_force = PB.get_if_available(vars.shelfarea_force, 1.0)
    @inbounds for i in cellrange.indices
        if pars.shelfareanorm[i] > 0.0
            ccrate_total = (shelfarea_force * pars.carbsedshallow[] * pars.shelfareanorm[i] 
                            * max(vars.OmegaAR[i]-1.0, 0.0)^1.7)
            ccrate = @PB.isotope_totaldelta(CIsotopeType, ccrate_total, vars.DIC_delta[i]+vars.D_mccb_DIC[])
            vars.shelf_Ccarb[i]            = ccrate  # NB: attribute :initialize_to_zero=true so we can skip cells with shelfareanorm=0

            fluxOceanfloorSolute.DIC[i]    -= ccrate
            fluxOceanfloorSolute.TAlk[i]   -= 2.0*PB.get_total(ccrate)
            fluxOceanBurial.Ccarb[i]       += ccrate
          
            PB.add_if_available(fluxOceanfloorSolute.Ca, i, -PB.get_total(ccrate))
        end
    end

    return nothing
end


"""
    ReactionBurialEffCarb

Deep ocean carbonate burial as a carbonate-saturation-state-dependent fraction of carbonate flux.

Parameterisation from [Caldeira1993](@cite), intended for use in ocean box models with a single deep ocean box.

Fraction of input carbonate flux buried `flys` is a function of oceanfloor carbonate saturation state.

Carbonate burial flux `flux_Ccarb` = `particulateflux_Ccarb` * `flys`, where Fraction of input carbonate flux 
buried `flys` is a function of carbonate saturation state.

Can be used with a spatially resolved model to provide a saturation-state dependent switch that allows burial of
oceanfloor carbonate flux only above the lysocline.

Parameter `burial_eff_function` sets the functional form used for burial fraction `flys`:
- "Caldeira1993": flys = 0.5*(1.0+tanh(k0([CO3] -k1)))
 (after [Caldeira1993](@cite), intended for use with the single deep ocean box in a 3-box ocean model)
- "OmegaCA": flys = 1 - 0.5 * erfc(m0*(Ω_CA - m1))  (burial efficiency a function of oceanfloor saturation state)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionBurialEffCarb{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParString("burial_eff_function", "Caldeira1993", allowed_values=["Caldeira1993", "OmegaCA"],
            description="functional form for burial efficiency"),
        PB.ParDouble("k0", 26.0, units="m3/mol",
            description="Caldeira1993: burial frac 'steepness' with CO3 concentration"),
        PB.ParDouble("k1", 0.11, units="mol/m3",
            description="Caldeira1993: CO3 concentration at burial frac 0.5"),
        PB.ParDouble("m0", 10.0, units="",
            description="OmegaCA: burial frac 'steepness' with OmegaCA"),
        PB.ParDouble("m1", 1.0, units="",
            description="OmegaCA: OmegaCA at burial frac 0.5"),

        PB.ParBoolVec("hascarbseddeep", Bool[],
            description="per box flag to enable (length=Domain size)"),

        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
    )


end

function PB.register_methods!(rj::ReactionBurialEffCarb)
  
    CIsotopeType = rj.pars.CIsotope[]
    PB.setfrozen!(rj.pars.CIsotope)
    @info "register_methods! $(PB.fullname(rj)) CIsotopeType=$CIsotopeType"
    
    fluxOceanfloorSolute = PB.Fluxes.FluxContrib(
        "fluxOceanfloor.soluteflux_", ["DIC::$CIsotopeType", "TAlk", "(Ca)"],
        alloptional=false,
    )

    fluxOceanBurial = PB.Fluxes.FluxContrib(
        "fluxOceanBurial.flux_", ["Ccarb::$CIsotopeType"],
    )

    vars = [
        PB.VarTarget("particulateflux_Ccarb", "mol yr-1", "input carbonate particulate flux",
            attributes=(:field_data=>CIsotopeType,)),        

        PB.VarProp("deep_Ccarb", "mol yr-1", "deep unresolved Ccarb burial",
            attributes=(:field_data=>CIsotopeType, :calc_total=>true,)),
        PB.VarProp("flys", "", "fraction of Ccarb export buried"),
    ]
    if CIsotopeType <: PB.AbstractIsotopeScalar
        push!(vars,
            PB.VarDep("ocean.oceanfloor.DIC_delta",          "per mil",  "d13C DIC"),
            PB.VarDepScalar("ocean.D_mccb_DIC",   "per mil",  
                "d13C marine calcite burial relative to ocean DIC"),
        )
    end

    function burial_eff_caldeira1993(pars, vars, i)
        return 0.5*(1.0 + tanh(pars.k0[]*(vars.CO3_conc[i] - pars.k1[])))
    end

    function burial_eff_omegaCA(pars, vars, i)
        return 1.0 - 0.5*erfc(pars.m0[]*(vars.OmegaCA[i] - pars.m1[]))
    end
    
    if rj.pars.burial_eff_function[] == "Caldeira1993"
        # ScalarData as we only want total, not isotopes (if any)
        push!(vars, PB.VarDep("ocean.oceanfloor.CO3_conc", "mol m-3", "CO3 concentration"; attributes=(:field_data=>PB.ScalarData,)))
        burial_eff_fn = burial_eff_caldeira1993
    elseif rj.pars.burial_eff_function[] == "OmegaCA"
        push!(vars, PB.VarDep("ocean.oceanfloor.OmegaCA", "", "calcite saturation"))
        burial_eff_fn = burial_eff_omegaCA
    else
        error("unknown burial_eff_function $(rj.pars.burial_eff_function[])")
    end    
    PB.setfrozen!(rj.pars.burial_eff_function)

    PB.add_method_do!(
        rj,
        do_burial_eff_carb,
        (
            PB.VarList_namedtuple_fields(fluxOceanfloorSolute),
            PB.VarList_namedtuple_fields(fluxOceanBurial),
            PB.VarList_namedtuple(vars),
        ),
        p = (CIsotopeType, burial_eff_fn), 
    )

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)

end

function PB.check_configuration(rj::ReactionBurialEffCarb, model::PB.Model)
    configok = true
 
    length(rj.pars.hascarbseddeep) == PB.get_length(rj.domain) || 
        (configok = false; @error "config error: length(carbseddeep) != Domain size")
 
    return configok
end

function do_burial_eff_carb(
    m::PB.ReactionMethod,
    pars,
    (fluxOceanfloorSolute, fluxOceanBurial, vars), 
    cellrange::PB.AbstractCellRange,
    deltat
)
    (CIsotopeType, burial_eff_fn) = m.p

    @inbounds for i in cellrange.indices 
        vars.flys[i] = burial_eff_fn(pars, vars, i)     
        if pars.hascarbseddeep[i]          
            ccrate = vars.flys[i] * vars.particulateflux_Ccarb[i]
        else
            ccrate = 0.0*vars.flys[i]*vars.particulateflux_Ccarb[i]
        end
        vars.deep_Ccarb[i]                  = ccrate

        # return all CaCO3 as solute Ca, CO3 except that which is buried
        ccremin = vars.particulateflux_Ccarb[i] - ccrate
        fluxOceanfloorSolute.DIC[i]    += ccremin
        fluxOceanfloorSolute.TAlk[i]   += 2.0*PB.get_total(ccremin)
        PB.add_if_available(fluxOceanfloorSolute.Ca, i, PB.get_total(ccremin))

        fluxOceanBurial.Ccarb[i]       += ccrate       
    end

    return nothing
end

# define a template instance so we know the type
const LININTERP_TEMPLATE = Interpolations.LinearInterpolation(
    [0.0, 1.0], 
    [0.0, 1.0],
    extrapolation_bc=Interpolations.Flat(),
)

"""
    ReactionBurialEffCorgP

Burial efficiency (fraction of `particulateflux` input) for Corg and P.

Input organic matter flux is given by `particulateflux_Corg, N, P`.  A fraction of `Corg` and `P` is buried
to output flux `fluxOceanBurial.flux_, Corg, P, Porg, PFe, Pauth`, where `P` is the total P burial flux and
`Porg, PFe, Pauth` are the three P burial mineral phases. The remainder of the input flux is
transferred to `reminflux_Corg, N, P`, where it would usually be linked to a `ReactionRemin` to be remineralized.
      
The fraction of `Corg` buried is given by a burial efficiency function, with options set by
Parameter `burial_eff_function`:
- `Prescribed`: Corg burial in cell `i` = Corg particulateflux * `BECorgNorm`*`Parameter BECorg[i]`
- `Ozaki2011`:  Corg burial in cell `i` =  Corg particulateflux * `BECorgNorm`*`burialEffCorg_Ozaki2011(sedimentation_rate)`
- `Dunne2007`: corg burial in cell `i` = Corg particulateflux * `BECorgNorm`*`burialEffCorg_Dunne2007(Corg particulateflux, Afloor)`
- `ConstantBurialRate`: Corg burial in cell `i` = `BECorgNorm`*`Parameter BECorg[i]` (independent of Corg flux)

It is also possible to set Parameter `FixedCorgBurialTotal`, in which case the ocean total Corg burial rate is fixed,
and the per-cell Corg burial fluxes calculated as above are then normalized to reach this global total rate.

P burial efficiency is defined as P:Corg ratios for components Porg, PFe, Pauth by parameters `BPorgCorg`, `BPFeCorg`, `BPauthCorg`.
If these are Vectors of length > 1, they define interpolated functions of oceanfloor `[O2]` on a grid defined by parameter `BPO2`.
If they are all Vectors of length 1, they define fixed Corg:P ratios.

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionBurialEffCorgP{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("BECorgNorm", 1.0, units="",
            description="overall normalization factor for Corg burial (or total Corg burial for ConstantBurialRate)"),
        PB.ParString("burial_eff_function", "Prescribed",
            allowed_values=["Prescribed", "Ozaki2011", "Dunne2007", "ConstantBurialRate"],
            description="Corg burial efficiency parameterisation (or ConstantBurialRate)"),
        PB.ParDoubleVec("BECorg", Float64[], units="",
            description="prescribed fraction seafloor Corg flux buried (or per-cell fraction of total Corg for ConstantBurialRate)"),
        PB.ParDouble("FixedCorgBurialTotal", NaN, units="mol C yr-1",
            description="if != NaN, fix total ocean Corg burial rate by renormalizing per-cell fluxes"),

        PB.ParDoubleVec("BPO2", [NaN], units="mol O2 m-3",
            description="[O2] points for interpolated oxygen-dependent P:Corg (length 1 for O2-independent P:Corg)"),
        PB.ParDoubleVec("BPorgCorg", [0.0], units="mol P (mol Corg)-1",
            description="P:Corg for organic P burial fraction at each [O2] (Vector length 1 for O2-independent P:Corg)"),
        PB.ParDoubleVec("BPFeCorg", [0.0], units="mol P (mol Corg)-1",
            description="P:Corg for Fe-associated P burial fraction at each [O2] (Vector length 1 for O2-independent P:Corg)"),
        PB.ParDoubleVec("BPauthCorg", [0.0], units="mol P (mol Corg)-1",
            description="P:Corg for CFA-associated P burial fraction at each [O2] (Vector length 1 for O2-independent P:Corg)"),

        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
    )

    fPorgCorg::typeof(LININTERP_TEMPLATE)   =  LININTERP_TEMPLATE
    fPFeCorg::typeof(LININTERP_TEMPLATE)    =  LININTERP_TEMPLATE
    fPauthCorg::typeof(LININTERP_TEMPLATE)  =  LININTERP_TEMPLATE
end

function PB.register_methods!(rj::ReactionBurialEffCorgP)

    CIsotopeType = rj.pars.CIsotope[]
    PB.setfrozen!(rj.pars.CIsotope)
    @info "register_methods! $(PB.fullname(rj)) CIsotopeType=$CIsotopeType"

    particulateflux = PB.Fluxes.FluxTarget(
        "particulateflux_", ["Corg::$CIsotopeType", "P", "N"],
        description="input particulate flux",
    )

    reminflux = PB.Fluxes.FluxContrib(
        "reminflux_", ["Corg::$CIsotopeType", "P", "N"],
        description="output unburied particulate flux",
        alloptional=false,
    )

    fluxOceanBurial = PB.Fluxes.FluxContrib(
        "fluxOceanBurial.flux_", ["Corg::$CIsotopeType", "Porg", "PFe", "Pauth", "P"],
    )

    vars = PB.VariableReaction[        
        PB.VarProp("burial_eff_Corg",      "", "Corg burial efficiency"),
        PB.VarDep("(sedimentation_rate)",      "m yr-1", "sedimentation rate"),
        PB.VarDep("(ocean.oceanfloor.O2_conc)", "mol m-3", "O2 concentration"),
        PB.VarDepStateIndep("(oceanfloor.Afloor)", "m^2",    "horizontal area of seafloor at base of box"),
    ]

    # Corg burial efficiency functions
    burial_eff_prescribed(pars, input_flux, vars, i) = pars.BECorgNorm[]*pars.BECorg[i]
                                                                                            
    burial_eff_Ozaki2011(pars, input_flux, vars, i) = pars.BECorgNorm[]*burialEffCorg_Ozaki2011(vars.sedimentation_rate[i])

    burial_eff_Dunne2007(pars, input_flux, vars, i) = pars.BECorgNorm[]*burialEffCorg_Dunne2007(input_flux.Corg[i], vars.Afloor[i])

    function burial_eff_forcedconstant(pars, input_flux, vars, i)
        forceconstCorg = pars.BECorgNorm[]*rj.pars.BECorg[i] # mol Corg yr-1 constant forced burial rate
        burial_eff = min(1.0, forceconstCorg/(PB.get_total(input_flux.Corg[i])+eps()))
        return burial_eff
    end

    if rj.pars.burial_eff_function[] == "Prescribed"
        burial_eff_fn = burial_eff_prescribed
    elseif rj.pars.burial_eff_function[] == "Ozaki2011"
        burial_eff_fn = burial_eff_Ozaki2011   
    elseif rj.pars.burial_eff_function[] == "Dunne2007"
        burial_eff_fn = burial_eff_Dunne2007        
    elseif rj.pars.burial_eff_function[] == "ConstantBurialRate"
        burial_eff_fn = burial_eff_forcedconstant
    else
        error("register_methods! $(PB.typename(rj)) $(PB.fullname(rj)) config error: "*
            "unknown burial_eff_function $(rj.pars.burial_eff_function[])")
    end
    PB.setfrozen!(rj.pars.burial_eff_function)

    PB.add_method_setup!(rj, setup_burial_eff_CorgP, (),)

    PB.add_method_do!(
        rj,
        do_burial_eff_CorgP,
        (
            PB.VarList_namedtuple_fields(particulateflux),
            PB.VarList_namedtuple_fields(reminflux),
            PB.VarList_namedtuple_fields(fluxOceanBurial),
            PB.VarList_namedtuple(vars)
        ),
        p = burial_eff_fn,
    )

    # PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end

function setup_burial_eff_CorgP(
    m::PB.ReactionMethod,
    pars,
    _,
    cellrange::PB.AbstractCellRange, 
    attribute_name,
)
    attribute_name == :setup || return

    rj = m.reaction

    if pars.burial_eff_function[] == "Prescribed"
        length(pars.BECorg) == PB.get_length(rj.domain) || 
            error("setup_burial_eff_CorgP! $(PB.fullname(rj)) config error: length(BECorg) != Domain size")
    end

    # Corg:P burial ratio functions
    all(length(pars.BPO2) .== length.((pars.BPorgCorg, pars.BPauthCorg, pars.BPFeCorg))) ||
        error("setup_burial_eff_CorgP! $(PB.fullname(rj)) config error: "*
            "lengths BPO2, BPorgCorg, BPauthCorg, BPFeCorg differ")
       
    if length(pars.BPO2) <= 1
        # oxygen-independent - duplicate values so that interpolation -> constant
        O2_vals = [-1.0, 1.0]
        BPorgCorg_vals = pars.BPorgCorg[[1, 1]]
        BPFeCorg_vals = pars.BPFeCorg[[1, 1]]
        BPauthCorg_vals = pars.BPauthCorg[[1, 1]]
    else
        O2_vals = pars.BPO2.v
        BPorgCorg_vals = pars.BPorgCorg.v
        BPFeCorg_vals = pars.BPFeCorg.v
        BPauthCorg_vals = pars.BPauthCorg.v
    end
    
    rj.fPorgCorg = Interpolations.LinearInterpolation(
        O2_vals, BPorgCorg_vals, 
        extrapolation_bc=Interpolations.Flat(),
    )
    rj.fPFeCorg = Interpolations.LinearInterpolation(
        O2_vals, BPFeCorg_vals, 
        extrapolation_bc=Interpolations.Flat(),
    )
    rj.fPauthCorg = Interpolations.LinearInterpolation(
        O2_vals, BPauthCorg_vals, 
        extrapolation_bc=Interpolations.Flat(),
    )

    return nothing
end

function do_burial_eff_CorgP(
    m::PB.ReactionMethod,
    pars,
    (particulateflux, reminflux, fluxOceanBurial, vars),
    cellrange::PB.AbstractCellRange, 
    deltat
)
    rj = m.reaction
    burial_eff_fn = m.p

    if length(pars.BPO2) > 1 
        !isnothing(vars.O2_conc) || error("oxygen-dependent C:P burial but no O2_conc variable")
    end

    # Calculate Corg burial efficiency first, to allow for optional normalization to pars.FixedCorgBurialTotal
    # zero of appropriate type for type stability
    burial_Corg_tot = zero(PB.get_total(particulateflux.Corg[1])*burial_eff_fn(pars, particulateflux, vars, 1))
    @inbounds for i in cellrange.indices 
        # Corg burial from burial efficiency
        vars.burial_eff_Corg[i] = burial_eff_fn(pars, particulateflux, vars, i)
        burial_Corg_tot += vars.burial_eff_Corg[i]*PB.get_total(particulateflux.Corg[i])
    end
    # optionally renormalize burial rate to fixed total Corg burial
    totnormfac = if isnan(pars.FixedCorgBurialTotal[])
        # 1.0 of appropriate type for type stability
        one(pars.FixedCorgBurialTotal[]/burial_Corg_tot)
    else
        length(cellrange.indices) == PB.get_length(rj.domain) ||
            error("FixedCorgBurialTotal can't be used with tiled cellrange")
        pars.FixedCorgBurialTotal[]/burial_Corg_tot
    end

    # Calculate burial fluxes
    @inbounds for i in cellrange.indices 
        burial_Corg = totnormfac*vars.burial_eff_Corg[i]*particulateflux.Corg[i]
        
        fluxOceanBurial.Corg[i] += burial_Corg
        # @Infiltrator.infiltrate
        reminflux.Corg[i] +=    particulateflux.Corg[i] - burial_Corg
        
        # no N burial
        reminflux.N[i]   += particulateflux.N[i]

        # P burial from burial efficiency
        O2_conc = PB.get_if_available(vars.O2_conc, i, -1.0) # not needed if oxygen-independent
        # @info "O2_conc = $(O2_conc), i = $(i)"
        BPorgCorg, BPFeCorg, BPauthCorg = rj.fPorgCorg(O2_conc), rj.fPFeCorg(O2_conc), rj.fPauthCorg(O2_conc)
        
        burial_Porg, burial_PFe, burial_Pauth = (BPorgCorg, BPFeCorg, BPauthCorg).*PB.get_total(burial_Corg)        
        burial_P        = burial_Porg + burial_PFe + burial_Pauth

        fluxOceanBurial.Porg[i]    += burial_Porg
        fluxOceanBurial.PFe[i]     += burial_PFe
        fluxOceanBurial.Pauth[i]   += burial_Pauth
        fluxOceanBurial.P[i]       += burial_P

        reminflux.P[i] += particulateflux.P[i] - burial_P
        if particulateflux.P[i] - burial_P < 0
            @warn "ReactionBurialEffCorgP burial_P exceeds particulateflux_P"
        end

    end

    return nothing
end

"""
    burialEffCorg_Ozaki2011(sedimentation_rate) -> BECorg

Organic carbon burial efficiency `BECorg` (fraction buried), [Ozaki2011](@cite)
as a function of `sedimentation_rate` (m/yr)
```jldoctest; setup = :(import PALEOocean)
julia> round(PALEOocean.Oceanfloor.Burial.burialEffCorg_Ozaki2011(2.369e-3), sigdigits=5) # 100m
0.26767
julia> round(PALEOocean.Oceanfloor.Burial.burialEffCorg_Ozaki2011(3.4152e-4), sigdigits=5) # 1000m
0.12335
```
"""
function burialEffCorg_Ozaki2011(sedimentation_rate)    
    # eqn (2) burial efficiency
              # SR in cm/yr
    BECorg = (100.0*sedimentation_rate)^0.4/2.1

    return BECorg
end

"""
    burialEffCorg_Dunne2007 -> BECorg

    %POC: percent organic carbon of surface sediments, 
    use log linear reegression of the data compilation by [Gedges and Keil1995](@cite)

    eqn(2) g cm^-2 a-1   g cm^-2 a-1
    %POC = FPOC_burial / Fmass_accumulation_rate * 100 = 12.0 * Fmass_accumulation_rate ^ 0.4

    To convert unit between sediment accumulation rates in velocity (cm/ka) to flux (g/cm a),
    the sediment porosity (φ = 0.7), density (ρ = 2.7) are used.

    BE = FPOC_burial/FPOC_bottom, the numerator is the POC burial flux, the denominator is flux of OC reeaching the sediments

    BE = 0.013 + 0.53 * (FPOC_bottom/(7.0+FPOC_bottom)) ^2 # [Dunne2007](@cite) and [Romaniello and Derry](@cite)
"""
function burialEffCorg_Dunne2007(
    Input_POC,  # mol/yr, could be a isotopelinear
    Afloor,     # m^2
    )
    # eqn (3) burial efficiency
    # the unit of FPOC_bottom (mol/yr) should be convert to (mmol C m^-2 d^-1)
    FPOC_bottom = PB.get_total(Input_POC) * 1000 / 365.25 / Afloor

    BECorg = 0.013 + 0.53 * (FPOC_bottom/(7.0+FPOC_bottom))^2
    
    return BECorg
end


end
