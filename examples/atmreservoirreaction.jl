"""
    temporary solution to make ReactionReservoirAtm available to PALEOocean examples

TODO move to PALEOboxes and remove from PALEOdev
"""
module AtmReservoirTODO

import PALEOboxes as PB
using PALEOboxes.DocStrings

"""
    ReactionReservoirAtm

A single scalar atmosphere reservoir (state variable).

Creates a state Variable `R` (units mol, with attribute `vfunction=VF_StateExplicit`) and 
(if parameter `const=false`) corresponding source-sink flux `R_sms` (units mol yr-1, with attribute `vfunction=VF_Deriv`).

Two associated Variables are also created:
- `pRatm`, partial pressure in units of bar or atm.  This is calculated by dividing `R` (mol) by
  parameter `moles1atm` (which has a default value for global modern Earth atmosphere).
- `pRPAL`, value normalized eg to a present-day modern Earth value. Calculated by dividing `R` (mol) by attribute `R:norm_value`,
  which should be set in the .yaml file to a representative value.

The local name prefix `R` should then be renamed using `variable_links:` in the configuration file.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionReservoirAtm{P} <: PB.AbstractReaction
    base::PB.ReactionBase
    
    pars::P = PB.ParametersTuple(
        PB.ParType(PB.AbstractData, "field_data", PB.ScalarData,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable isotopes and specify isotope type"),
        PB.ParBool("const", false,
            description="true to provide constant value with no _sms Variable"),
        PB.ParDouble("moles1atm", PB.Constants.k_moles1atm, units="mol",
            description="moles for a pressure of 1 Earth atm.  Default value $(PB.Constants.k_moles1atm) for global modern Earth atmosphere")
    )

    norm_value::Float64 = NaN
end

function PB.register_methods!(rj::ReactionReservoirAtm)
 
    @info "register_methods! ReactionReservoirAtm $(PB.fullname(rj)) field_data=$(rj.pars.field_data[])"
    PB.setfrozen!(rj.pars.field_data)

    if rj.pars.const[]
        R = PB.VarPropScalar(      "R", "mol", "scalar constant reservoir", attributes=(:field_data=>rj.pars.field_data[],))
        R_sms = PB.VarTarget(     "R_sms", "mol yr-1", "constant reservoir source-sinks", attributes=(:field_data=>rj.pars.field_data[],))
    else
        R = PB.VarStateExplicitScalar("R", "mol", "scalar reservoir", attributes=(:field_data=>rj.pars.field_data[],))
        R_sms = PB.VarDerivScalar(     "R_sms", "mol yr-1", "scalar reservoir source-sinks", attributes=(:field_data=>rj.pars.field_data[],))
    end
    PB.setfrozen!(rj.pars.const)

    # sms variable not used by us, but must appear in a method to be linked and created
    PB.add_method_do_nothing!(rj, [R_sms])

    vars = [
        R,
        PB.VarPropScalar(         "pRnorm", "", "scalar atmosphere reservoir normalized eg to present-day value"),
        PB.VarPropScalar(         "pRatm", "atm", "scalar atmosphere reservoir partial pressure"),
    ]
       
    if rj.pars.field_data[] <: PB.AbstractIsotopeScalar
        push!(vars, 
            PB.VarPropScalar(     "R_delta", "per mil", "scalar atmosphere reservoir isotope delta"))
    end

    # callback function to store Variable norm during setup
    function setup_callback(m, attribute_value, v, vdata)
        v.localname == "R" || error("setup_callback unexpected Variable $(PB.fullname(v))")
        if attribute_value == :norm_value
            m.reaction.norm_value = PB.value_ad(PB.get_total(vdata[]))
        end
        return nothing
    end
    # set filterfn to force setup even if R is constant, not a state Variable
    PB.add_method_setup_initialvalue_vars_default!(rj, [R], filterfn = v->true, setup_callback=setup_callback)  

    PB.add_method_do!(
        rj,
        do_reservoir_atm,
        (PB.VarList_namedtuple(vars), ),
    )

    PB.add_method_initialize_zero_vars_default!(rj)
    
    return nothing
end

function do_reservoir_atm(m::PB.AbstractReactionMethod, pars, (vars, ), cr::PB.AbstractCellRange, deltat)
    rj = m.reaction

    vars.pRnorm[]  = PB.get_total(vars.R[])/rj.norm_value
    vars.pRatm[]  = PB.get_total(vars.R[])/pars.moles1atm[]

    if hasfield(typeof(vars), :R_delta)
        vars.R_delta[] = PB.get_delta(vars.R[])
    end

    return nothing
end

end # module