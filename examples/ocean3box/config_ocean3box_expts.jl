include("../atmreservoirreaction.jl") # temporary solution to make ReactionReservoirAtm available

"test cases and examples for 3 box ocean"
function config_ocean3box_expts(model, expts)

    ###########################
    # configure expt
    ############################

    for expt in expts        
        println("Add expt: ", expt)
        if expt == "baseline"
            # defaults
        elseif expt == "fastexchange"
            react_airseaCO2 = PB.get_reaction(model, "oceansurface", "airsea_CO2")
            PB.setvalue!(PB.get_parameter(react_airseaCO2, "piston"), 3.1e3)  # m/day 'fast exchange'
        elseif expt == "slowexchange"
            react_airseaCO2 = PB.get_reaction(model, "oceansurface", "airsea_CO2")
            PB.setvalue!(PB.get_parameter(react_airseaCO2, "piston"), 3.1e-3)  # m/day 'slow exchange'
        elseif expt == "killbio"
            # disable biology at t=0
            react_enable_bioprod = PB.get_reaction(model, "global", "force_enable_bioprod")
            PB.setvalue!(PB.get_parameter(react_enable_bioprod, "force_times"), [-1e30, 0.0, 1.0e-16, 1e30]) 
            PB.setvalue!(PB.get_parameter(react_enable_bioprod, "force_values"), [1.0, 1.0, 0.0, 0.0]) 
        elseif expt == "lowO2"
            # start with low oxygen to test marine sulphur system
            PB.set_variable_attribute!(model, "atm", "O2", :initial_value, 0.1*3.71e19)
        elseif expt == "lowSO4"
            # start with low SO4 to test methane cycling
            PB.set_variable_attribute!(model, "ocean", "SO4", :initial_value, 100e-3)  # ~100 uM       
        else
            error("unrecognized expt='$(expt)'")
        end
    end

    return nothing
end
