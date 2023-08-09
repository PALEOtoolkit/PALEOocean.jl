"Ocean production Reactions"
module BioProd

import PALEOboxes as PB
using PALEOboxes.DocStrings

import PALEOaqchem

# import Infiltrator # Julia debugger
    
"""
    ReactionBioProdPrest

P-limited biological production, configurable to restore P to specified level or to consume fraction of nutrients

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionBioProdPrest{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParIntVec("bioprod",  Int[], allowed_values=[0,1,2,3],
            description="production type (per cell): 0 - none, 1 - restore to absolute conc, 2 - consume fraction of nutrients supplied, 3 - restore to fraction of ocean mean"),
        PB.ParDoubleVec("bioprodval",  Float64[], units="m-3 or none", 
            description="conc or frac corresponding to 'bioprod'"),
        PB.ParDouble("rCorgPO4",  106.0, units="", 
            description="Corg:P Redfield ratio of organic matter produced"),
        PB.ParDouble("rNPO4",     16.0, units="",
            description="N:P Redfield ratio of organic matter produced"),
        PB.ParDouble("rCcarbCorg", 0.0, units="",
            description="ratio of Ccarb to Corg produced"),
        PB.ParDouble("trest",  0.1, units="yr",
            description="restoring timescale"),
        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
    )
end

function PB.register_methods!(rj::ReactionBioProdPrest)

    CIsotopeType = rj.pars.CIsotope[]
    PB.setfrozen!(rj.pars.CIsotope)

    vars = [
        PB.VarDep("volume",               "m^3",      "ocean cell volume"),
        PB.VarDepScalar("volume_total",   "m^3",      "ocean total volume"),
        PB.VarDep("P_conc",               "mol m^-3", "total P concentration"),
    
        PB.VarProp("%reaction%Prod_Corg", "mol yr-1", "organic carbon production rate",
                                    attributes=(:field_data=>CIsotopeType, :calc_total=>true,)),
        PB.VarProp("%reaction%Prod_Ccarb", "mol yr-1", "carbonate production rate",
                                    attributes=(:field_data=>CIsotopeType, :calc_total=>true,)),
        PB.VarContrib("P_sms",            "mol yr-1", "total dissolved P source minus sink"),
        PB.VarContrib("(O2_sms)",         "mol yr-1", "O2 source minus sink"),
        PB.VarContrib("(TAlk_sms)",       "mol yr-1", "TAlk source minus sink"),
        PB.VarContrib("(DIC_sms)",        "mol yr-1", "DIC source minus sink", 
                                    attributes=(:field_data=>CIsotopeType,)),
    
        PB.VarDepScalar("(global.enable_bioprod)", "", "optional forcing, =0.0 to disable, !=0.0 to enable"),
    
        PB.VarDepScalar("(global.PELCALC)", "", "optional forcing for pelagic calcification"),
    ]

    # add optional Variables needed for par_bioprod options
    if any(rj.pars.bioprod .== 2)
        push!(vars, PB.VarDep("P_transport_input",  "mol yr-1", "transport P input"))
    end        
    if any(rj.pars.bioprod .== 3)
        push!(vars, PB.VarDepScalar("P_total",      "mol",      "total P"))
    end
    PB.setfrozen!(rj.pars.bioprod)

    # add Variable required for isotopes
    if CIsotopeType <: PB.AbstractIsotopeScalar
        push!(vars,
            PB.VarDep("DIC_delta",          "per mil",  "d13C DIC"),
            PB.VarDepScalar("D_mccb_DIC",   "per mil",  "d13C marine calcite burial fractionation relative to ocean DIC"),
            PB.VarDepScalar("D_B_mccb_mocb","per mil",  "d13C fractionation between marine organic and calcite burial"),
        )
    end


    prod_vars = PB.Fluxes.FluxContrib(
        "prod_", ["P", "N", "Corg::$CIsotopeType", "Ccarb::$CIsotopeType"], 
    )
   
    PB.add_method_do!(
        rj,
        do_bio_prod_Prest,
        (PB.VarList_namedtuple(vars), PB.VarList_namedtuple_fields(prod_vars)),
        p = CIsotopeType
    )

    PB.add_method_do_totals_default!(rj)

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end

# Override to validate parameters
function PB.check_configuration(rj::ReactionBioProdPrest, model::PB.Model)
    configok = true
    
    length(rj.pars.bioprod) == PB.get_length(rj.domain) ||
        (configok = false; @error "$(PB.fullname(rj)) invalid Parameter 'bioprod' length != Domain size")
    
    length(rj.pars.bioprodval) == PB.get_length(rj.domain)   || 
        (configok = false; @error "$(PB.fullname(rj)) invalid Parameter 'bioprodval' length != Domain size")
    
    return configok
end

function do_bio_prod_Prest(m::PB.ReactionMethod, pars, (vars, prod_vars), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction
    CIsotopeType = m.p

    for i in cellrange.indices   
        ProdTotP = zero(vars.P_conc[i])    # use zero() to get correct type eg if using AD
        if isnothing(vars.enable_bioprod) || vars.enable_bioprod[] != 0
            if pars.bioprod[i] == 0
                # no production       
                # ProdTotP = 0.0     
            elseif pars.bioprod[i] == 1
                # restore to absolute P level
                # mol yr-1        mol m-3                                  m^3          / yr
                ProdTotP += max(vars.P_conc[i] - pars.bioprodval[i], 0.0)*vars.volume[i]/pars.trest[]
            elseif pars.bioprod[i] == 2
                # consume a fraction of nutrients supplied by circulation into this box
                ProdTotP += pars.bioprodval[i] * vars.P_transport_input[i] 
            elseif pars.bioprod[i] == 3
                # restore to fraction of ocean mean P level
                mean_P_conc = vars.P_total[]/vars.volume_total[]
                restore_P_conc = mean_P_conc * pars.bioprodval[i]                
                # mol yr-1        mol m-3                                  m^3          / yr
                ProdTotP += max(vars.P_conc[i] - restore_P_conc, 0.0)*vars.volume[i]/pars.trest[]
            else           
                error("bioprod=$(pars.bioprod[i]) not supported")
            end
        end
        ProdTotN        = ProdTotP * pars.rNPO4[]

        ProdTotCorg_total = ProdTotP * pars.rCorgPO4[]
        ProdTotCcarb_total = ProdTotCorg_total*pars.rCcarbCorg[]*PB.get_if_available(vars.PELCALC, 1.0)
        if CIsotopeType <: PB.AbstractIsotopeScalar        
            # calculate d13C
            delta13C_prod_ccarb = vars.DIC_delta[i] + vars.D_mccb_DIC[]
            delta13C_prod_corg = delta13C_prod_ccarb - vars.D_B_mccb_mocb[]                
        end
        ProdTotCorg         = @PB.isotope_totaldelta(CIsotopeType, ProdTotCorg_total, delta13C_prod_corg)
        vars.Prod_Corg[i]   = ProdTotCorg                 
        ProdTotCcarb        = @PB.isotope_totaldelta(CIsotopeType, ProdTotCcarb_total, delta13C_prod_ccarb)
        vars.Prod_Ccarb[i]  = ProdTotCcarb

        # corresponding O2, TAlk uptake
        (uptakeOrgO2, uptakeAlk) = PALEOaqchem.O2AlkUptakeRemin(ProdTotCorg_total, (ProdTotN, 0, 0),  ProdTotP,  ProdTotCcarb_total)

        # all state variable sms except P_sms are optional to allow use with eg a P only configuration
        vars.P_sms[i] += -ProdTotP
        PB.add_if_available(vars.O2_sms, i, -uptakeOrgO2)
        PB.add_if_available(vars.TAlk_sms, i, -uptakeAlk)
        PB.add_if_available(vars.DIC_sms, i, -(ProdTotCorg + ProdTotCcarb))

        PB.add_if_available(prod_vars.P, i, ProdTotP)
        PB.add_if_available(prod_vars.N, i, ProdTotN)
        PB.add_if_available(prod_vars.Corg, i, ProdTotCorg)
        PB.add_if_available(prod_vars.Ccarb, i, ProdTotCcarb)       

    end

    return nothing
end

"""
    ReactionBioProdMMPop

Ocean phytoplankton biological production.

Configurable to represent oxygenic photosynthesizers with P, N limitation, nitrogen fixers, anoxygenic photosynthesis limited
by electron donor availability.

Export production or production is represented as a combination of limiting factors:

`population_size x nutrient_limitation x light_limitation x temperature_limitation x electron_donor_limitation`

`population_size` may be either implicit (either a constant or ∝ nutrient concentration, generating GENIE-like parameterisations of 
export production), or represented explicitly as a state variable (in which case production is accumulated into a state variable `phytP_conc`,
and `k_grazeresprate` defines a background loss rate that is exported).

Export production is than partitioned into DOM flux (components `domprod_P`, `domprod_N`, `domprod_Corg`)
and particulate flux (components `partprod_P`, `partprod_N`, `partprod_Corg`, `partprod_Ccarb`) fractions according to parameter `k_nuDOM`.

See [Kriest2010](@cite) for a comparison of models of this type.

## Production functional forms
## `population_size`
Set by `k_poptype` parameter:
- `Constant`: `k_uPO4` (mol P / m-3 / yr) (represents export production)
- `Nutrient`: `k_mu*P_conc` (mol P / m-3 / yr) (represents export production)
- `Pop` : `k_mu*phytP_conc` (mol P / m-3 / yr) (represents growth of explicit population `phtyP_conc`)

## `nutrient_limitation` 
Set by `k_nuttype` parameter:
- `PO4MM`: `P_conc / (P_conc + k_KPO4)`
- `PO4NMM`: P, N limited phytoplankton export production: cf GENIE `2N2T_PO4MM_NO3`, [Fennel2005](@cite)
- `PO4NMMNfix`: P limited nitrogen-fixer export production: cf GENIE `2N2T_PO4MM_NO3`, [Fennel2005](@cite)

## `light_limitation`
Set by `k_lightlim` parameter:
- `fixed`: constant `k_Irel`
- `linear`: GENIE-like form `k_Irel*insol/PB.Constants.k_solar_presentday` where `insol` (W m-2) is provided insolation in each cell
- `MM`: MITgcm-like saturating form `zInsol/(k_Ic + zInsol)` where `zInsol = k_Irel*insol` and `insol` is provided PAR in each cell
- `QE`: Saturating light limitation of rate vs (local) PAR `k_Irel*insol`, derived from photosynthetic QE `k_alphaQE` and chl absorption cross section `k_thetaChlC`.
  See eg [Geider1987](@cite) for summary of notation and unit conversions.

## `temperature_limitation`
Set by `k_templim` parameter:
- `Constant`: constant value 1.0
- `Eppley`: Eppley curve, normalized to 1.0 at 0 deg C,  `exp(0.0633*(temp - PB.Constants.k_CtoK))`

## `electron_donor_limitation`
Set by `k_edonor` parameter:
- `H2O`: constant 1.0, no electron-donor limitation on production
- `H2S`: `H2S_conc / (H2S_conc + k_KH2S)` production limited by H2S concentration

# Examples

| `k_poptype` | `k_nuttype` | `k_lightlim` | `k_templim` | `k_edonor` | Reference                           |
|:------------|:------------|:-------------|:------------|:-----------|:------------------------------------|
| Constant    | PO4MM       | linear       | Constant    | H2O        | P and light limited export production, as used by GENIE [Ridgwell2007](@cite) |
| Pop         | PO4MM       | QE           | Eppley      | H2O        | P and light limited phytoplankton population  |

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionBioProdMMPop{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        # Parameters controlling stoichiometry of organic matter and carbonate produced
        PB.ParDouble("rCorgPO4",  106.0, units="", 
            description="Corg:P Redfield ratio of organic matter produced"),
        PB.ParDouble("rNPO4",     16.0, units="",
            description="N:P Redfield ratio of organic matter produced"),
        PB.ParDouble("rCcarbCorg", 0.0, units="",
            description="ratio of Ccarb to Corg produced"),
        PB.ParBool("rCcarbCorg_fixed", true,
            description = "Ccarb:Corg rain ratio true for fixed, false for sat. state dependent"),
        PB.ParDouble("k_r0",  0.044372,
            description = "initial rain-ratio for sat. state dependent rain ratio"),
        PB.ParDouble("k_eta", 0.8053406,
            description= "exponent for sat. state dependent rain ratio"),

        PB.ParDouble("nuDOM",  0.66, units="", 
            description="fraction of production to DOM reservoir"),
        PB.ParDouble("depthlimit", -200.0, units="m",
            description="depth limit for production"),

        # Parameters for population growth rate model   
        PB.ParString("k_poptype", "Constant", allowed_values=["Constant", "Nutrient", "Pop"],
            description = "population / growth rate model"),
        # Parameters controlling (max) growth rate and losses
        PB.ParDouble("k_uPO4", 2.0e-3, units="mol P / m-3 / yr",
            description="for k_poptype = 'Constant': max rate, constant ie of form k_O_uPO4 * (light, nut etc)"),
        PB.ParDouble("k_mu", NaN, units="1/yr",
            description="for k_poptype = 'Nutrient', 'Pop': max prod/growth rate (at 0C if templim=='Eppley')"),
        PB.ParDouble("k_grazeresprate", NaN, units="1/yr",
            description="for k_poptype = 'Pop' imposed const loss (mortality) rate"),
    
        # Parameters for temperature limitation
        PB.ParString("k_templim", "Constant", allowed_values=["Constant", "Eppley"],
            description="temperature limitation factor"),

        # Parameters for  light limitation
        PB.ParString("k_lightlim", "linear", allowed_values=["linear", "MM", "fixed", "QE"],
            description="Light limitation function"),
        PB.ParDouble("k_Irel",  1.0, 
            description="multiplier for forcing-supplied insolation"),
        PB.ParDouble("k_Ic", 30.0, units="W/m^2",
            description="saturating irradiance for 'MM' case"),
        PB.ParDouble("k_alphaQE", 7.0, units="mgC/mgChl/Wpar m^-2/d-1",
            description="chla-specific initial slope of the photosynthesis-light curve for lightlim='QE'"),
        PB.ParDouble("k_thetaChlC", 0.03, units="mg Chl / mgC",
            description="Chl:Corg ratio for explicit population k_poptype=Pop"),
        PB.ParDouble("k_epsilonChl", 0.012, units="m^2/mg Chl",
            description="chl absorption coeff (for self shielding) for k_poptype='Pop'"),
    
        # Parameters for nutrient limitation
        PB.ParString("k_nuttype", "PO4MM", allowed_values=["PO4MM", "PO4NMM", "PO4NMMNfix"],
            description="Nutrient limitation / nitrogen fixation function"),
        PB.ParDouble("k_KPO4", NaN,   units="mol P m-3 ",
            description="limitation at low [P] MM half-max constant"),
        PB.ParDouble("k_KN", 0.0, units="mol NO3+NH4 m-3",
            description="limitation at low nitrogen"),
        PB.ParDouble("k_prefNH3", 10.0, units="",
            description="preference for ammonia over nitrate"),
        
        # Parameters controlling electron donor
        PB.ParString("k_edonor", "H2O",     allowed_values=["H2O", "H2S"],
            description="electron donor (H2O for oxygenic phototroph"),
        PB.ParDouble("k_KH2S", 1.0e-3,   units="mol H2S m-3",
            description="limitation at low [H2S] MM half-max constant"),

        # Isotopes
        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
    )
end

function PB.register_methods!(rj::ReactionBioProdMMPop)
    
    CIsotopeType = rj.pars.CIsotope[]
    PB.setfrozen!(rj.pars.CIsotope)

    grid_vars = [
        PB.VarDep("volume",               "m^3",      "ocean cell volume"),  
        PB.VarDep("zupper",               "m",        "cell depth (-ve)"),
    ]

    vars = [        
        PB.VarDepColumn("oceansurface.open_area_fraction","","fraction of area open to atmosphere"),
       
        PB.VarDep("P_conc",               "mol m^-3", "total P concentration"),
        PB.VarDep("(OmegaCA)",            "",         "calcite saturation state"),

        PB.VarDepScalar("(global.rate_bioprod)", "", "optional forcing, multiplier for productivity"),
        PB.VarDepScalar("(global.PELCALC)", "", "optional forcing for pelagic calcification"),
    
                                    # set :initialize_to_zero as do_react may only cover some cells (within depth range)
        PB.VarProp("%reaction%Prod_Corg", "mol yr-1", "organic carbon production rate",
            attributes=(:field_data=>CIsotopeType, :initialize_to_zero=>true, :calc_total=>true)),
        PB.VarProp("%reaction%Prod_Ccarb", "mol yr-1", "carbonate production rate",
            attributes=(:field_data=>CIsotopeType, :initialize_to_zero=>true, :calc_total=>true)),
    
        PB.VarContrib("P_sms",            "mol yr-1", "total dissolved P source minus sink"),
        PB.VarContrib("(TAlk_sms)",       "mol yr-1", "TAlk source minus sink"),
        PB.VarContrib("(DIC_sms)",        "mol yr-1", "DIC source minus sink",
            attributes=(:field_data=>CIsotopeType,)),
    ]

    append!(vars,
         PB.Fluxes.FluxContrib("partprod_", ["P", "N", "Corg::$CIsotopeType", "Ccarb::$CIsotopeType"])
    )
    append!(vars,
        PB.Fluxes.FluxContrib("domprod_", ["P", "N", "Corg::$CIsotopeType", "Ccarb::$CIsotopeType"])
    )

    if CIsotopeType <: PB.AbstractIsotopeScalar
        push!(vars,
            PB.VarDep("DIC_delta",          "per mil",  "d13C DIC"),
            PB.VarDepScalar("D_mccb_DIC",   "per mil",  "d13C marine calcite burial relative to ocean DIC"),
            PB.VarDepScalar("D_B_mccb_mocb","per mil",  "d13C fractionation between marine organic and calcite burial"), 
        )
    end
    
    # get functions for edonor, nutrient, light limitation
    # these add any required VariableReactions to vars
    (f_edonor_lim, edonor_ratestoich) = get_edonor_functions(vars, rj.pars.k_edonor[])
    PB.setfrozen!(rj.pars.k_edonor)
    PB.add_method_do!(rj, edonor_ratestoich)
    
    f_nuttype                       = get_nuttype_function(rj.pars.k_nuttype[])
    PB.setfrozen!(rj.pars.k_nuttype)

    (f_rate_poptype, f_loss_poptype)= get_rate_poptype_functions(rj, grid_vars, vars, rj.pars.k_poptype[])
    PB.setfrozen!(rj.pars.k_poptype)

    f_lightlim                      = get_lightlim_function(vars, rj.pars.k_lightlim[])
    PB.setfrozen!(rj.pars.k_lightlim)

    f_templim                       = get_templim_function(vars, rj.pars.k_templim[])
    PB.setfrozen!(rj.pars.k_templim)
    
    PB.add_method_do!(
        rj,
        do_bio_prod_MM_pop,
        (
            PB.VarList_namedtuple(grid_vars),
            PB.VarList_namedtuple(vars),
        ),
        p = (CIsotopeType, f_edonor_lim, f_nuttype, f_rate_poptype, f_loss_poptype, f_lightlim, f_templim),
    )

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end


function do_bio_prod_MM_pop(
    m::PB.ReactionMethod,
    pars,
    (grid_vars, vars),
    cellrange::PB.AbstractCellRange,
    deltat
)

    CIsotopeType, f_edonor_lim, f_nuttype, f_rate_poptype, f_loss_poptype, f_lightlim, f_templim = m.p

    for (isurf, colindices) in cellrange.columns
        for i in colindices
            if grid_vars.zupper[i] > pars.depthlimit[]
                # dimensionless (0-1)
                (lim_nut, (fracTNH3, fracNO3, fracNfix)) = f_nuttype(pars, vars, i)

                # dimensionless (0-1)
                lim_edonor = f_edonor_lim(pars, vars, i)
              
                # mol P/m^3/yr                          
                rate_pop  = f_rate_poptype(pars, vars, i)
             
                # dimensionless
                lim_temp = f_templim(pars, vars, i)

                # dimensionless
                lim_Ifac = f_lightlim(pars, vars, i, pars.k_mu[]*lim_temp)

                # optional forcing for production rate
                rate_bioprod          = PB.get_if_available(vars.rate_bioprod, 1.0)

                # mol P yr-1
                ProdTotP = (
                    rate_bioprod 
                    *rate_pop  # mol P m-3 yr-1
                    *lim_temp
                    *lim_nut
                    *lim_edonor
                    *lim_Ifac
                    *vars.open_area_fraction[isurf]
                    *grid_vars.volume[i] # m^3
                )

                # mol P yr-1
                LossTotP = f_loss_poptype(pars, vars, i, ProdTotP)        

                # stoichiometry
                NtoP        = pars.rNPO4[]
                CorgtoP     = pars.rCorgPO4[]
                if pars.rCcarbCorg_fixed[]
                    CcarbtoCorgfac    = pars.rCcarbCorg[]
                else
                    CcarbtoCorgfac    = pars.k_r0[]*max(vars.OmegaCA[i]-1, 0.0)^pars.k_eta[]
                end
                CcarbtoP    = CcarbtoCorgfac * pars.rCorgPO4[]*PB.get_if_available(vars.PELCALC, 1.0)
            
                # corresponding O2, TAlk uptake for Corg only (to subtract from tracer _sms)
                # -ve
                # @Infiltrator.infiltrate
                (uptakeOrgO2eqtoP, uptakeAlkOrgtoP) = PALEOaqchem.O2AlkUptakeRemin(
                                                    CorgtoP, 
                                                    NtoP.*(fracNO3, fracTNH3, fracNfix),  
                                                    1.0,  
                                                    0.0)
            
                if CIsotopeType <: PB.AbstractIsotopeScalar
                    # calculate d13C
                    delta13C_ccarb = vars.DIC_delta[i] + vars.D_mccb_DIC[]
                    delta13C_corg = delta13C_ccarb - vars.D_B_mccb_mocb[]                
                end
                ProdTotCorg         = @PB.isotope_totaldelta(CIsotopeType, ProdTotP*CorgtoP, delta13C_corg)
                vars.Prod_Corg[i]   = ProdTotCorg
                LossTotCorg         = @PB.isotope_totaldelta(CIsotopeType, LossTotP*CorgtoP, delta13C_corg)
                ProdTotCcarb        = @PB.isotope_totaldelta(CIsotopeType, ProdTotP*CcarbtoP, delta13C_ccarb)
                vars.Prod_Ccarb[i]  = ProdTotCcarb
                LossTotCcarb        = @PB.isotope_totaldelta(CIsotopeType, LossTotP*CcarbtoP, delta13C_ccarb)

                # Losses -> particulates (export)
                PB.add_if_available(vars.partprod_P,    i, (1-pars.nuDOM[])*LossTotP)
                PB.add_if_available(vars.partprod_N,    i, (1-pars.nuDOM[])*LossTotP*NtoP)
                PB.add_if_available(vars.partprod_Corg, i, (1-pars.nuDOM[])*LossTotCorg)
                PB.add_if_available(vars.partprod_Ccarb,i, (1-pars.nuDOM[])*LossTotCcarb)       

                # Losses -> DOM
                PB.add_if_available(vars.domprod_P,    i, pars.nuDOM[]*LossTotP)
                PB.add_if_available(vars.domprod_N,    i, pars.nuDOM[]*LossTotP*NtoP)
                PB.add_if_available(vars.domprod_Corg, i, pars.nuDOM[]*LossTotCorg)
                # NB: no Ccarb -> DOM, so add this back to DIC/TAlk immediately as a 'short circuit'
                CcarbtoDOM = pars.nuDOM[]*LossTotCcarb

                # Tendencies (nutrients etc consumed)
                # all state variable sms except P_sms are optional to allow use with eg a P only configuration
                vars.P_sms[i] += -ProdTotP
                vars.edonorO2eq[i] = -uptakeOrgO2eqtoP*ProdTotP
                    
                PB.add_if_available(vars.TAlk_sms, i, -uptakeAlkOrgtoP*ProdTotP - 2.0*PB.get_total(ProdTotCcarb) + 2.0*PB.get_total(CcarbtoDOM))
                PB.add_if_available(vars.DIC_sms, i, -(ProdTotCorg + ProdTotCcarb) + CcarbtoDOM)
                # TODO nitrogen

            end
        end            
    end

    return nothing
end

"""
Michaelis-Menton like resource limitation
    lim_MM = P_conc/(P_conc+k_KPO4)
"""
@inline function lim_MM(val, halfmax)
    lim = max(PB.get_total(val), 0.0)/(max(PB.get_total(val), 0.0) + halfmax)
    return lim
end

function get_nuttype_function(nuttype::AbstractString)

    "smooth function approximating step function:
    at large k,  h(x<0) = 0, h(x>0) = 1"
   hsmooth(x, k=100.0) = 1.0/(1.0+exp(-k*x))

   # Nutrient limitation functions:  nut_X(rj, vars, ixd) -> (nutlim, (fracTNH3, fracNO3, fracNfix))
   "P limitation only"    
   function nut_PO4MM(pars, vars, i)
       lim_Pfac = lim_MM(vars.P_conc[i], pars.k_KPO4[])
       return (lim_Pfac, (0.0, 1.0, 0.0))
   end

   "P,N limited phytoplankton: GENIE 2N2T_PO4MM_NO3 (cf Fennel etal 2005)"
   function nut_PO4NMM(pars, vars, i)
        lim_Pfac = lim_MM(vars.P_conc[i], pars.k_KPO4[])
        
        cTNH3 = max(vars.TNH3_conc[i], 0.0)
        cNO3 = max(vars.NO3_conc[i], 0.0)
        lim_Nfac = (cTNH3 + cNO3)/(cTNH3 + cNO3 + pars.k_KN[])
        fracTNH3 = pars.k_prefNH3[]/(pars.k_O_prefNH3[]*cTNH3 + cNO3 + eps())
        fracNO3 = 1.0-fracTNH3
        fracNfix = 0.0

       return (min(lim_Pfac, lim_Nfac), (fracTNH3, fracNO3, fracNfix))
   end

   "P limited nitrogen fixer: GENIE 2N2T_PO4MM_NO3 (cf Fennel etal 2005)"
   function nut_PO4NMMNfix(pars, vars, i)
        lim_Pfac = lim_MM(vars.P_conc[i], pars.k_KPO4[])
       
        # Nfix where NO3 + TNH3 < abs(obj.k_O_KN)
        cN = max(vars.TNH3[i], 0.0) + max(vars.NO3_conc[i], 0.0)
        Nlowfac = hsmooth((pars.k_O_KN[] - cN)/pars.k_KN[])  # > 0 when N fix allowed
       
        # Nfix where N:P < redfield
        Nredfac = hsmooth(16.0 - cN./max(vars.P_conc[i],eps()))
                  
        # Combine both criteria: ->1 where N fix permitted
        lim_Nfac = Nlowfac*Nredfac

        return (lim_Pfac*lim_Nfac, (0.0, 0.0, 1.0))
   end

   if nuttype == "PO4MM"
        return nut_PO4MM
   elseif nuttype == "PO4NMM"
        return nut_PO4NMM
   elseif nuttype == "PO4NMMNfix"
        return nut_PO4NMMNfix
   else
        error("unknown nuttype ", nuttype)
   end
end

"Stoichiometry for H2O as e- donor"
const default_edonorH2O = PB.RateStoich(
    PB.VarProp("edonorO2eq", "mol O2eq yr-1", "O2eq e- donor consumption (H2O) by oxygenic photosynthesis",
        # set :initialize_to_zero as do_react may only cover some cells (within depth range)
        attributes=(:initialize_to_zero=>true, :calc_total=>true)),
    ((1.0, "O2"),),
    sms_prefix="", sms_suffix="_sms",
    processname="production",
)

"Stoichiometry for H2S as e- donor"
const default_edonorH2S = PB.RateStoich(
    PB.VarProp("edonorO2eq", "mol O2 eq yr-1", "O2eq e- donor consumption (H2S) by anoxygenic photosynthesis",
            # set :initialize_to_zero as do_react may only cover some cells (within depth range)
            attributes=(:initialize_to_zero=>true, :calc_total=>true)),
    ((-0.5, "H2S::Isotope"), (0.5, "SO4::Isotope"), (-1.0, "TAlk")),
    deltavarname_eta = ("H2S_delta", -35.0),   # TODO constant fractionation 35 per mil ??
    sms_prefix="", sms_suffix="_sms",
    processname="production",
)   

function get_edonor_functions(vars, edonor::AbstractString)

    function edonorlim_nolimit(pars, vars, i)
        return 1.0  # no edonor limitation on production
    end

    function edonorlim_H2S(pars, vars, i)
        lim_H2Sfac = lim_MM(vars.H2S_conc[i], pars.k_KH2S[])
        return lim_H2Sfac
    end

    if edonor == "H2O"
        edonor_ratestoich = deepcopy(default_edonorH2O)
        push!(vars, edonor_ratestoich.ratevartemplate) # edonorO2eq
        return (edonorlim_nolimit, edonor_ratestoich)
    elseif edonor == "H2S"
        edonor_ratestoich = deepcopy(default_edonorH2S)
        push!(vars, edonor_ratestoich.ratevartemplate) # edonorO2eq
        # ScalarData as we only want concentration, not isotopes (if any)
        push!(vars, PB.VarDep("H2S_conc", "mol m^-3", "total H2S concentration", attributes=(:field_data=>PB.ScalarData,)))
        return (edonorlim_H2S, edonor_ratestoich)
    else
        error("unknown edonor ", edonor)
    end    

    return nothing
end

function get_rate_poptype_functions(rj, grid_vars, vars, poptype::AbstractString)

    # (mol P/m^3/yr, 1/yr)

    function rate_constant(pars, vars, i)
        return pars.k_uPO4[]
    end

    function rate_nutrient(pars, vars, i)
        return pars.k_mu[]*max(vars.P_conc[i], 0.0)
    end

    "explicit phytoplankton population"
    function rate_pop(pars, vars, i)
        # 1/yr max growth rate       
        rate_pop2 = pars.k_mu[]*vars.phytP_conc[i]   # mol P / m^3 / yr 
        # TODO type inference fails (generating spurious allocations and slow code)
        # if function name same as (return?) variable name !?
        return rate_pop2  # see above - use a different variable name to work around a Julia 1.6 bug
    end

    # mol P yr-1
    function poploss_nopop(pars, vars, i, ProdTotP)
        return ProdTotP
    end

    # mol P yr-1
    function poploss_const(pars, vars, i, ProdTotP)
        # mol P yr-1   mol P    *  yr-1
        LossTotP = vars.phytP[i]*pars.k_grazeresprate[]
        vars.phytP_sms[i] += ProdTotP - LossTotP
        return LossTotP
    end

    if poptype == "Constant"
        return (rate_constant, poploss_nopop)
    elseif poptype == "Nutrient"
        return (rate_nutrient, poploss_nopop)
    elseif poptype == "Pop"
        phytP = PB.VarStateExplicit("%reaction%phytP", "mol P", "phytoplankton population",
            attributes=(:calc_total=>true,)
        )
        phytP_sms = PB.VarDeriv("%reaction%phytP_sms", "mol P yr-1", "phytoplankton population source - sink")
        phytP_conc = PB.VarProp("%reaction%phytP_conc", "mol P m-3", "phytoplankton population concentration",
            attributes=(:advect=>true, :vertical_movement=>0.0, :specific_light_extinction=>0.0)
        )
        pop_vars = [
            phytP,
            phytP_sms,
            phytP_conc,
            PB.VarContrib("(opacity)",        "m-1",  "total optical opacity"),                           
        ]
        PB.add_method_setup_initialvalue_vars_default!(rj, [phytP])
        PB.add_method_do!(rj, do_bio_prod_MM_popconc, (PB.VarList_namedtuple(grid_vars), PB.VarList_namedtuple(pop_vars)) )

        push!(vars, PB.VarContrib(phytP_sms))
        push!(vars, PB.VarDep(phytP))
        push!(vars, PB.VarDep(phytP_conc))

        return (rate_pop, poploss_const)
    else
        error("unknown poptype ", poptype)
    end
end

function get_templim_function(vars, templimtype::AbstractString)
    function templim_const(pars, vars, i)
        return 1.0
    end

    function templim_eppley(pars, vars, i)
        # 1.0 at 0 deg C, Eppley curve
        return exp(0.0633*(vars.temp[i] - PB.Constants.k_CtoK))
    end

    if templimtype == "Constant"
        return templim_const
    elseif templimtype == "Eppley"
        push!(vars, PB.VarDep("temp",                "K",    "temperature"))
        return templim_eppley
    else
        error("unknown templim_type ", templim_type)
    end
end

function get_lightlim_function(vars, lightlim::AbstractString)
    "GENIE - like proportional to insolation"
    function lightlim_linear(pars, vars, i, muMaxT)
        lim_Ifac = pars.k_Irel[]*vars.insol[i]/PB.Constants.k_solar_presentday
        return lim_Ifac
    end

    "MITgcm - like saturating form"
    function lightlim_MM(pars, vars, i, muMaxT)
        zInsol = pars.k_Irel[]*vars.insol[i]  # PAR at depth z including wc attenuation
        lim_Ifac = zInsol/(pars.k_Ic[] + zInsol)
        return lim_Ifac
    end

    "specified value"
    function lightlim_fixed(pars, vars, i, muMaxT)
        lim_Ifac = pars.k_Irel[]
        return lim_Ifac
    end

    """
        Saturating light limitation of rate vs (local) irradiance,
        derived from photosynthetic QE and chl absorption cross section.
        See eg Geider (1987) New Phytol. for summary of notation and unit conversions.
            phiQE = 0.09;  // mol C (mol photon PAR)-1
            EINSTEINPERWM2 = 5.0; // conversion factor umol photons m-2 s-1 per W PAR m-2
            gC (g Chla)-1 (Wm-2)-1 d-1
            alphaQE = 1.62e-5*EINSTEINPERWM2*SECPERDAY*(phiQE/0.09) = 7.0 gC (g Chla)-1 (Wm-2)-1 d-1
        
                            thetaChlC  * par                  
        check value: alphaQE*0.03       * 4.8 [*(epsilonChl/0.015) ?]  = 1.0 d-1
    """
    function lightlim_QE(pars, vars, i, muMaxT)
        par = pars.k_Irel[]*vars.insol[i]   # W/m^2 PAR at depth z 
        #                    mgC/mgChl/Wpar m^-2/d         d / yr            * mgChl/mgC * Wpar m^2 / (1/yr)  
        QEfac_lightlim = pars.k_alphaQE[]*PB.Constants.k_daypyr*pars.k_thetaChlC[]*par/muMaxT
        lim_Ifac = (1 - exp(-QEfac_lightlim))
        return lim_Ifac
    end

    var_insol = PB.VarDep("insol",                "W m-2",    "photosynthetic radiative flux")
    if lightlim == "linear"
        push!(vars, var_insol)
        return lightlim_linear
    elseif lightlim == "MM"
        push!(vars, var_insol)
        return lightlim_MM
    elseif lightlim == "fixed"
        return lightlim_fixed
    elseif lightlim == "QE"
        push!(vars, var_insol)
        return lightlim_QE
    else
        error("unknown lightlim $lightlim")
    end

end



function do_bio_prod_MM_popconc(
    m::PB.ReactionMethod,
    pars,
    (grid_vars, pop_vars),
    cellrange::PB.AbstractCellRange,
    deltat
)

    @inbounds for i in cellrange.indices
        pop_vars.phytP_conc[i] = pop_vars.phytP[i]/grid_vars.volume[i]            
        # mg Chl m-3       mg Chl / mgC        mol C / mol P          mol P m-3        mg C / mol C    
        phytmgchlm3 = pars.k_thetaChlC[]*pars.rCorgPO4[]*pop_vars.phytP_conc[i]*12000
        # m-1                    m^2/mg Chl        mg Chl / m^3                 
        PB.add_if_available(pop_vars.opacity, i, pars.k_epsilonChl[] * phytmgchlm3)
    end
  
    return nothing
end

end
