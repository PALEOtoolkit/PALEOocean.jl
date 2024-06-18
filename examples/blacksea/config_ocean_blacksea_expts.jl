function config_ocean_blacksea_expts(
    model, expts
)
    for expt in expts
        if expt == "baseline"
            # baseline configuration

        elseif length(expt) == 5 && expt[1] == "setpar"
            # generic parameter set (setpar, <domain>, <reaction>, <parname>, <parvalue)
            _, domname, reactname, parname, parvalue = expt
            rct = PB.get_reaction(model, domname, reactname)
            !isnothing(rct) || error("expt $expt Reaction $domname.$reactname not found")
            par = PB.get_parameter(rct, parname)
            PB.setvalue!(par, parvalue)
            
        
        else
            error("unknown expt ", expt)
        end
    end

    return model
end
