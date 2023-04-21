module AirSeaExchange

import PALEOocean

import PALEOboxes as PB
using PALEOboxes.DocStrings

# import Infiltrator
                       
const ScCoeff = (
    #       A       -B          +C         -D
    CO2 = (2073.1, -125.62,    +3.6276,    -0.043219),  # Wanninkhof (1992) Table A1
    CH4 = (2039.2, -120.31,    +3.4209,    -0.040437),  # Wanninkhof (1992) Table A1
    O2  = (1638.0, -81.53,     +1.483,     -0.008004),  # Sarmiento & Gruber (2006) Table 3.3.1, which cites Keeling etal (1998)
)

"""
    Sc = SchmidtNumber(temperatureK, coeff)

Seawater (salinity 35 per mil) Schmidt number (dimensionless)
julia> round(PALEOocean.Oceansurface.AirSeaExchange.ScO2(273.15+15), sigdigits=14)
721.7115
"""
function SchmidtNumber(temperatureK, coeff)
    
    ScTminC = 0.0   # temperature limits for Schmidt factors 
    ScTmaxC = 30.0

    T_degC = temperatureK - PB.Constants.k_CtoK
    T_degC = max(T_degC, ScTminC)
    T_degC = min(T_degC, ScTmaxC)
    
    # Sc = A - B T + C T^2 - D T^3

    Sc = coeff[1] + coeff[2]*T_degC + coeff[3]*T_degC^2 + coeff[4]*T_degC^3

    return Sc
end

ScCO2(temperatureK) = SchmidtNumber(temperatureK, ScCoeff.CO2)
ScCH4(temperatureK) = SchmidtNumber(temperatureK, ScCoeff.CH4)

"""
    ScO2(temperatureK)

Seawater (salinity 35 per mille) Schmidt number (dimensionless)

Check values:
```jldoctest; setup = :(import PALEOocean)
julia> round(PALEOocean.Oceansurface.AirSeaExchange.ScO2(273.15+15), sigdigits=14)
721.7115
```
"""
ScO2(temperatureK)  = SchmidtNumber(temperatureK, ScCoeff.O2)

"""
    satH2O(temperatureK, salinity, applylimit=true)

Saturation H2O vapour pressure (atm)
Sarmiento & Gruber Panel 3.2.1, citing Weiss and Price (1980)

Check values:
```jldoctest; setup = :(import PALEOocean)
julia> round(PALEOocean.Oceansurface.AirSeaExchange.satH2O(273.15+15, 35), sigdigits=15)
0.0164958943614262
```
 """
function satH2O(temperatureK, salinity, applylimit=true)
    
    # default behaviour is to apply limits, can be disabled for testing
    if applylimit
        # TODO a guess at limits
        temperatureK = max(temperatureK, PB.Constants.k_CtoK - 2.0)
        temperatureK = min(temperatureK, PB.Constants.k_CtoK + 50.0)
        
        salinity = max(salinity,0.0)
        salinity = min(salinity,40.0)
    end
    
    satH2O = exp(24.4543 - 67.4509*(100.0/temperatureK) - 4.8489*log(temperatureK/100.0)
                 - 0.000544*salinity)

    return satH2O
end

"""
    kwWann1992(Sc, Uwind) -> kw

Gas transfer velocity Wanninkhof (1992)

# Arguments:
 - `Sc: Schmidt number
 - `Uwind`: (m/s) wind velocity at 10m above sea level

# Return
 - `kw`: (m/day) piston velocity

# Check values:
```jldoctest; setup = :(import PALEOocean)
julia> round(PALEOocean.Oceansurface.AirSeaExchange.kwWann1992(720.0, 10.0), sigdigits=15)
7.12325768170716
```
"""
function kwWann1992(Sc, Uwind)
    
    k_O_a = 0.31   # scalar: Wanninkhof (1992)
    
    piston = (0.01*             # cm -> m
              24*               # hr -> day
              k_O_a*            # scalar
              Uwind^2*          # windspeed (m s-1)
              (Sc/660.0)^-0.5)  # Schmidt number

    return piston
end


                        
const SolubilityCoeff = (
    #       A1            A2          A3          B1          B2          B3               
    # beta from fit in Weiss (1970) (also Wanninkhof (1992) Table A2)
    O2 = (-58.3877,     85.8079,    23.8439,    -0.034892,  0.015568,   -0.0019387),
    # Wanninkhof (1992) Table A2
    CH4= (-68.8862,     101.4596,   28.7314,    -0.076146,  0.043970,   -0.0068672),
)

"""
    gas_solubility(coeff, temperatureK, salinity, applylimit=true) -> sol

solubility (mol / l / P_moist atm) ie function of (moist) air partial pressure

T and salinity limits are for O2, discussion in text Weiss (1970)
dataset (Fig 3) includes 1 <= TdegC <= 36 
Table 7 includes -1 <= TdegC <= 40, 0 <= sal <= 40
discussion in text describes calculation of isotherms for -2 <= T <= 40
"""
function gas_solubility(coeff, temperatureK, salinity, applylimit=true)
    # default behaviour is to apply limits, can be disabled for testing
    if applylimit
        temperatureK = max(temperatureK, PB.Constants.k_CtoK - 2.0)
        temperatureK = min(temperatureK, PB.Constants.k_CtoK + 40.0);
        
        salinity = max(salinity, 0.0)
        salinity = min(salinity, 40.0)
    end

    # beta  units [volume of gas at STP / unit volume of solution / atm gas] 
    beta = exp(
                coeff[1] +
                coeff[2]*(100.0/temperatureK) +
                coeff[3]*log(temperatureK/100.0) +
                salinity*(
                    coeff[4] +
                    coeff[5]*(temperatureK/100.0) +
                    coeff[6]*((temperatureK/100.0)^2))
                )
            
    # convert to volumetric solubility assuming ideal gas behaviour
    sol = beta/PB.Constants.k_molVolIdealGas   # mol gas / litre solution / atm partial pressure gas

    return sol
end

"""
    solO2Wann92(temperatureK, salinity, applylimit=true)

solubility (mol / l / P_moist atm) ie function of (moist) air partial pressure.

Check values:
```jldoctest; setup = :(import PALEOocean)
julia> round(PALEOocean.Oceansurface.AirSeaExchange.solO2Wann92(273.15+25, 30), sigdigits=14)
0.0010692031999441
julia> round(PALEOocean.Oceansurface.AirSeaExchange.solO2Wann92(273.15+25, 30), sigdigits=14)
0.0010692031999441
```
"""
solO2Wann92(temperatureK, salinity, applylimit=true) =
    gas_solubility(SolubilityCoeff.O2, temperatureK, salinity, applylimit)

solCH4Wann92(temperatureK, salinity, applylimit=true) =
    gas_solubility(SolubilityCoeff.CH4, temperatureK, salinity, applylimit)

"""
    frac_airsea_CO2(tempK) -> (eps_k, eps_eqb) 

Kinetic and equilibrium fractionation delta13C for air-sea transfer
"""
function frac_airsea_CO2(tempK)
    # Fractionation from Zhang, Quay, Wilbur (1995)
    eps_k                      = -0.9;  # kinetic fractionation
    eps_eqb                    = 10.53-0.105*(tempK - PB.Constants.k_CtoK)  # equilibrium fractionation

    return (eps_k, eps_eqb)
end

"""
    frac_airsea_none(tempK) -> (0.0, 0.0) 

No kinetic or equilibrium fractionation for air-sea transfer
"""
frac_airsea_none(tempK) = (0.0, 0.0)


"""
    ReactionAirSea(Sc_function, sol_function, frac_function, gas_name)

Not called directly - use as `ReactionAirSeaO2`, `ReactionAirSeaCO2`, `ReactionAirSeaCH4`, `ReactionAirSeaFixedSolubility`

Atmosphere to ocean gas X flux, mol/yr.
"""
Base.@kwdef mutable struct ReactionAirSea{P} <: PB.AbstractReaction
    base::PB.ReactionBase
        
    pars::P = PB.ParametersTuple(
        PB.ParBool("solubility_fixed", false, 
            description="use fixed solubility"),
        PB.ParDouble("sol_fix_henry_coeff", NaN, units="mol l-1 atm-1",
            description="Henry's law coefficient for solubility_fixed=true"),
        PB.ParBool("moistair",    true,
            description="apply correction for moist air"),
        PB.ParBool("piston_fixed",true,
            description="use fixed piston velocity"),
        PB.ParDouble("piston",    NaN, units="m d-1",
            description="fixed piston velocity"),
        PB.ParDouble("TempKmin",  -Inf, units="K",
            description="GENIE bug compatibility - lower limit on temperature for solubility"),
    )
    
    gas_name::String

    "function returning gas X Schmidt number  Sc_function(tempK) -> Sc"
    Sc_function
    "function returning gas X solubility sol_function(tempK, sal) -> sol (mol/l/atm)
        = nothing to derive solubility from oceansurface X_conc and pXatm"
    sol_function
    "function returning gas X isotope kinetic and eqb fractionation frac_function(tempK) -> (eps_k, eps_eqb) (per mil)"
    frac_function = nothing

    isotope_name::String    = ""
end
   
PB.create_reaction(::Type{ReactionAirSea}, base::PB.ReactionBase) =
    error("ReactionAirSea cannot be created directly - use a gas-specific version eg ReactionAirSeaO2")

"""
    ReactionAirSeaO2

See [`ReactionAirSea`](@ref)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
abstract type ReactionAirSeaO2 <: PB.AbstractReaction end
PB.create_reaction(::Type{ReactionAirSeaO2}, base::PB.ReactionBase) =
    create_ReactionAirSea(base, ScO2, solO2Wann92, false, nothing, "O2", "")

"""
    ReactionAirSeaCO2

See [`ReactionAirSea`](@ref)

# Parameters
$(PARS)
"""
abstract type ReactionAirSeaCO2 <: PB.AbstractReaction end
PB.create_reaction(::Type{ReactionAirSeaCO2}, base::PB.ReactionBase) =
    create_ReactionAirSea(base, ScCO2, nothing, false, frac_airsea_CO2, "CO2", "CIsotope")

"""
    ReactionAirSeaCH4

See [`ReactionAirSea`](@ref)

# Parameters
$(PARS)
"""
abstract type ReactionAirSeaCH4 <: PB.AbstractReaction end
PB.create_reaction(::Type{ReactionAirSeaCH4}, base::PB.ReactionBase) =
    create_ReactionAirSea(base, ScCH4, solCH4Wann92, false, frac_airsea_none, "CH4", "CIsotope")

"""
    ReactionAirSeaFixedSolubility

See [`ReactionAirSea`](@ref)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
abstract type ReactionAirSeaFixedSolubility <: PB.AbstractReaction end
PB.create_reaction(::Type{ReactionAirSeaFixedSolubility}, base::PB.ReactionBase) =
    create_ReactionAirSea(base, nothing, nothing, true, nothing, "X", "")


"Create new instance of reaction, parameterised by functions for Schmidt number, solubility, and fractionation.
Set `sol_fixed_value` to `true` to create generic Reaction with fixed solubility"
function create_ReactionAirSea(
    base, 
    Sc_function, 
    sol_function,
    solubility_fixed::Bool,
    frac_function, 
    gas_name, 
    isotope_name
)
    rj = ReactionAirSea(
        base=base,
        Sc_function=Sc_function,
        sol_function=sol_function,
        frac_function=frac_function,
        gas_name=gas_name,
        isotope_name=isotope_name
    )

    if solubility_fixed
        PB.setvalueanddefault!(rj.pars.solubility_fixed, true, freeze=true)
    end

    return rj
end


function PB.register_methods!(rj::ReactionAirSea)
    
    if isempty(rj.isotope_name)
        IsotopeType = PB.ScalarData
    else
        _, IsotopeType = PB.split_nameisotope("::"*rj.isotope_name, rj.external_parameters)
        @info "$(PB.fullname(rj)) IsotopeType=$IsotopeType from isotope_name=$(rj.isotope_name) in external_parameters"
    end
    
    vars = [
        PB.VarDep("Asurf",                              "m^2",      "horizontal area of oceansurface"),
        PB.VarDep("open_area_fraction",                 "",         "fracton of surface open to atmosphere (0-1.0)"),
       
        PB.VarDepScalar("pXatm"=>"atm.p"*rj.gas_name*"atm",   "atm",      "gas X atmospheric partial pressure",
            attributes=(:field_data=>PB.ScalarData,),),  # just the total, not isotopic composition (if any)    
        PB.VarDep("X_conc"=>"ocean.oceansurface."*rj.gas_name*"_conc", 
                                                        "mol m-3",  "ocean concentration [gas X]",
            attributes=(:field_data=>PB.ScalarData,),),  # just the total, not isotopic composition (if any)
        
        PB.VarContrib("flux_X"=>"fluxAtmtoOceansurface.flux_"*rj.gas_name,
                                                        "mol yr-1", "air -> sea gas X flux",
            attributes=(:field_data=>IsotopeType,),),
    ]

    if !rj.pars.piston_fixed[]
        push!(vars, PB.VarDep("wind_speed",             "m s-1", "Surface wind speed at 10 m"))
    else
        PB.setfrozen!(rj.pars.piston_fixed)
    end

    if !rj.pars.solubility_fixed[]
        push!(vars,
            PB.VarDep("ocean.oceansurface.temp",            "K",        "ocean surface temperature"),
            PB.VarDep("ocean.oceansurface.sal",             "psu",      "ocean salinity"),
        )
        if isnothing(rj.sol_function)
            push!(vars,
                PB.VarDep("pXocean"=>"ocean.oceansurface.p"*rj.gas_name, "atm", "ocean gas X partial pressure")
            )
        end
    end
        
    if IsotopeType <: PB.AbstractIsotopeScalar
        push!(vars,
            PB.VarDepScalar("Xatm_delta"=>"atm."*rj.gas_name*"_delta", "per mil",  "gas X atmosphere isotope delta"),
            PB.VarDep("Xocean_delta"=>"ocean.oceansurface."*rj.gas_name*"_delta", "per mil",  "gas X isotope delta"),
        )
    end

   
    PB.add_method_do!(
        rj, 
        do_air_sea_flux,
        (PB.VarList_namedtuple(vars), ),
        # add Sc_function, sol_function so these are supplied as arguments in order to maintain type stability
        # also add IsotopeType so this can be dispatched on
        p = (IsotopeType, rj.Sc_function, rj.sol_function, rj.frac_function)
    )

    return nothing
end

# Override to validate parameters
function PB.check_configuration(rj::ReactionAirSea, model::PB.Model)
    configok = true

    if rj.pars.piston_fixed[]
        rj.pars.piston[] >= 0 || 
            (configok = false; @error "$(PB.fullname(rj)) invalid Parameter 'piston'=$(rj.pars.piston[])")
    else
        !isnothing(rj.Sc_function) ||
            (configok = false; @warn " Parameter 'piston_fixed'==false but no Schmidt function supplied")   
    end
        
   
    return configok
end

function do_air_sea_flux(
    m::PB.ReactionMethod,
    pars,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    (IsotopeType, Sc_function, sol_function, frac_function) = m.p

    @inbounds for i in cellrange.indices            
        # piston velocity 
        if pars.piston_fixed[]
            # fixed piston velocity
            piston = pars.piston[] # m/day
        else
            # Schmidt number
            Sc = Sc_function(vars.temp[i])                        
            # calculate piston velocity (m/day) (see Wanninkhof 1992)              
            piston = kwWann1992(Sc, vars.wind_speed[i]) # m/day
        end
        # convert to volume exch/year, with consideration of sea-ice cover
        pfac = piston*PB.Constants.k_daypyr*vars.Asurf[i]*vars.open_area_fraction[i] # m^3 yr^-1
                                                    
        # Calculate solubility (mol/kg-sw/atm) for surface boxes                  
        if pars.solubility_fixed[]
            solconstX = pars.sol_fix_henry_coeff[]
        else
            tempKeff = max(vars.temp[i], pars.TempKmin[])  # set optional floor on temperature for GENIE bug compatibility
            if isnothing(sol_function)
                # mol/l/atm     mol m-3  m3 l-1 / atm
                solconst = PB.get_total(vars.X_conc[i])*1e-3/vars.pXocean[i]  
            else
                solconst = sol_function(tempKeff, vars.sal[i])  # mol/l/atm pp of X
            end
            if pars.moistair[]
                # correct for moist air with pO2atm interpreted as mixing ratio of O2 in dry air
                satH2Oval = satH2O(tempKeff, vars.sal[i])  # saturation H2O vapour pressure              
                solconstX = solconst*(1-satH2Oval) # mol / l / dry air mixing ratio 
            else
                # omit correction - ie wrong
                solconstX = solconst
            end
        end

        if !(IsotopeType <: PB.AbstractIsotopeScalar)
            # mol/yr          = m^3/yr * mol/l/atm* l/m^3 * atm                  mol/m^3)
            vars.flux_X[i]  += pfac*(solconstX*1000.0*vars.pXatm[] - max(vars.X_conc[i],0))
        else
            # isotopes - need to explicitly calculate flux in both directions to account for kinetic fractionation
            eps_k, eps_eqb = frac_function(vars.temp[i])

            fluxAtoOcean = PB.isotope_totaldelta(
                IsotopeType,
                pfac*solconstX*1000.0*vars.pXatm[], 
                vars.Xatm_delta[]+eps_eqb+eps_k
            )
            fluxOceantoA = PB.isotope_totaldelta(
                IsotopeType,
                pfac*PB.get_total(vars.X_conc[i]), 
                vars.Xocean_delta[i]+eps_k
            )
            vars.flux_X[i] += fluxAtoOcean - fluxOceantoA
        end
    end

    return nothing
end


end
