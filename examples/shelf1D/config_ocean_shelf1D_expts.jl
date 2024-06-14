
 function config_ocean_shelf1D_expts(
    model, expts,
)

    for expt in expts        
        println("Add expt: ", expt)
        if expt == "baseline"
            # defaults
        elseif length(expt) == 2 && expt[1] == "setO2"
            # ("setO2", <pO2 PAL> start with low oxygen to test marine sulphur system
            PB.set_variable_attribute!(model, "atm", "O2", :initial_value, expt[2]*3.71e19)
        elseif length(expt) == 2 && expt[1] == "setSO4"
            # ("setSO4", <SO4 conc mol m-3>) start with low SO4 to test methane cycling           
            PB.set_variable_attribute!(model, "ocean", "SO4", :initial_value, expt[2])
        elseif length(expt) == 2 && expt[1] == "setH2S"
            # ("setH2S", <H2S conc mol m-3>) start with high H2S
            PB.set_variable_attribute!(model, "ocean", "H2S", :initial_value, expt[2])
        elseif length(expt) == 2 && expt[1] == "setCH4"
            # ("setCH4", <CH4 conc mol m-3>) start with high CH4
            PB.set_variable_attribute!(model, "ocean", "CH4", :initial_value, expt[2])
        elseif expt == "QE"
            react_bioprod = PB.get_reaction(model, "ocean", "bioprod")
                         
            k_lightlim = react_bioprod.pars.k_lightlim[]
            k_lightlim == "QE" ||
                error("ocean.bioprod Reaction k_lightlim $k_lightlim != QE (change parameter in yaml file)")

            PB.setvalue!(react_bioprod.pars.k_Irel,  0.4) # multiplier to correct to PAR            
            PB.setvalue!(react_bioprod.pars.k_alphaQE,  7.0) # gC/gChl/Wpar m^-2/d-1
            PB.setvalue!(react_bioprod.pars.k_thetaChlC, 0.01) # TODO require a low value otherwise continuous growth through winter

        elseif expt == "nooxphot"
            react_bioprod = PB.get_reaction(model, "ocean", "bioprod")                                   
            PB.setvalue!(react_bioprod.pars.k_Irel,  0.0) # no PAR -> disable oxygenic photosynthesis
         
        else
            error("unrecognized expt='$(expt)'")
        end
    end

    return model
end

