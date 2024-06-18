
function config_ocean_romglb_expts(
        model, expts; 
)
    for expt in expts
        if expt == "baseline"
            # baseline configuration

        elseif expt == "set_carbonate_config_79box"
            shelfareanorm = zeros(79)
            shelfareanorm[14] = 1.0 # low lat (:gyre) surface box only
            r_shelfcarb = PB.get_reaction(model, "oceanfloor", "shelfcarb")
            PB.setvalue!(r_shelfcarb.pars.shelfareanorm, shelfareanorm)
    
            hascarbseddeep = ones(Bool, 79)
            hascarbseddeep[[1, 14, 47]] .= false # disable deposition in surface boxes ('shelves')
            r_deepcarb = PB.get_reaction(model, "oceanfloor", "deepcarb")
            PB.setvalue!(r_deepcarb.pars.hascarbseddeep, hascarbseddeep)
        
        else
            error("unknown expt ", expt)
        end
    end

    return model
end
