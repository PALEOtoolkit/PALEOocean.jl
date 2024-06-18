module SedimentationRate_dev


import PALEOboxes as PB
using PALEOboxes.DocStrings

using SpecialFunctions
using Interpolations

import Infiltrator



"""
    sedRate_Tromp1995(depth) -> sedimentation_rate

`sedimentation_rate` (m/yr) at `depth` (m, -ve), [Tromp1995](@cite)
eg used by [Ozaki2011](@cite) 

# Examples:
Check value at depth 1000m
 ```jldoctest
julia> round(Main.SedimentationRate_dev.sedRate_Tromp1995(-1000.0), sigdigits=5)
0.00034152
julia> round(Main.SedimentationRate_dev.Burial.sedRate_Tromp1995(-100.0), sigdigits=5)
0.002369
```
"""
function sedRate_Tromp1995(depth)    
    # eqn (3) sed rate: (NB: log = log10)

    sedimentation_rate = NaN*one(depth)

    if (depth < -5400)
        sedimentation_rate = zero(depth)
    elseif depth < 0
        sedimentation_rate = 1e-2*10^(erfcinv(-depth/2700)-2.1)
    end
  
    return sedimentation_rate
end

"""
    ReactionSedimentationRate_dev

Sedimentation rate parameterized from water depth

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_SETUP)
"""
Base.@kwdef mutable struct ReactionSedimentationRate_dev{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParString("sed_rate_function", "Tromp1995", allowed_values=["Tromp1995"],
            description="sedimentation rate parameterisation"),
    )

end

function PB.register_methods!(rj::ReactionSedimentationRate_dev)
    
    sr_Tromp1995(pars, vars, i) = sedRate_Tromp1995(vars.zlower[i])

    if rj.pars.sed_rate_function[] == "Tromp1995"
        sr_func = sr_Tromp1995
    else
        error("configuration error: unknown sed_rate_function $(rj.pars.sed_rate_function[])")
    end
    PB.setfrozen!(rj.pars.sed_rate_function)

    vars = [
        PB.VarDepStateIndep("ocean.oceanfloor.zlower", "m",    "depth of lower surface of box (m)"),
        PB.VarPropStateIndep("sedimentation_rate",      "m yr-1", "sedimentation rate"),
    ]
    
    PB.add_method_setup!(
        rj,
        setup_sedimentation_rate,
        (PB.VarList_namedtuple(vars), ),
        p = sr_func,
    )
    return nothing
end

function setup_sedimentation_rate(
    m::PB.ReactionMethod,
    pars,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    attribute_name
)
    attribute_name == :setup || return nothing

    sr_func = m.p

    @inbounds for i in cellrange.indices         
        vars.sedimentation_rate[i]  = sr_func(pars, vars, i)
    end

    return nothing
end


end
