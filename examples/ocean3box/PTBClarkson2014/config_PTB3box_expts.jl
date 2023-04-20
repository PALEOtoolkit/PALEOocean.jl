
"test cases and examples for 3 box ocean"
function config_PTB3box_expts(baseconfig, expts)

    if  baseconfig=="oaopencarb"
        # Open atmosphere-ocean with silicate carbonate weathering input and carbonate burial

        model = PB.create_model_from_config(
            joinpath(@__DIR__, "PALEO_examples_ocean3box_cfg.yaml"), "ocean3box_oaopencarb_base"
        )
            
    elseif  baseconfig=="Co2HOmLWCpp"
        # Clarkson (2014) CO2Hi

        model = PB.create_model_from_config(
            joinpath(@__DIR__, "PALEO_examples_PTB3box_cfg.yaml"), "ocean3box_oaopencarb_base"
        )    

        # constant_degassing = PB.get_reaction(model, "atm", "constant_degassing")
        # PB.setvalue!(constant_degassing.pars.perturb_totals, 11.8e12.*[1.0, 1.0])  # TODO d13C

        land_Bergman2004 = PB.get_reaction(model, "land", "land_Bergman2004")
        PB.setvalue!(land_Bergman2004.pars.k_silw, 2.40e12)
        PB.setvalue!(land_Bergman2004.pars.k17_oxidw, 5.00e12)   
        
        B.set_variable_attribute!(model, "ocean", "B_conc", :initial_delta,  36.8)

        shelfcarb = PB.get_reaction(model, "oceanfloor", "shelfcarb")
        PB.setvalue!(shelfcarb.pars.carbsedshallow, 18.43e12)

    elseif  baseconfig=="Co2LOmHWC4pp"
        # Clarkson (2014)  CO2Lo

        model = PB.create_model_from_config(
            joinpath(@__DIR__, "PALEO_examples_PTB3box_cfg.yaml"), "ocean3box_oaopencarb_base")
            

        # constant_degassing = PB.get_reaction(model, "atm", "constant_degassing")
        # PB.setvalue!(constant_degassing.pars.perturb_totals, 11.8e12.*[1.0, 1.0])  # TODO d13C

        land_Bergman2004 = PB.get_reaction(model, "land", "land_Bergman2004")
        PB.setvalue!(land_Bergman2004.pars.k_silw, 6.60e12)
        PB.setvalue!(land_Bergman2004.pars.k17_oxidw, 5.92e12)       

        PB.set_variable_attribute!(model, "ocean", "B_conc", :initial_delta, 34.0)

        shelfcarb = PB.get_reaction(model, "oceanfloor", "shelfcarb")
        PB.setvalue!(shelfcarb.pars.carbsedshallow, 1.44e12)

    else
        error("unrecognized baseconfig='$(baseconfig)'")
    end

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
        elseif expt == "killbioEP2"
            # disable ocean biology at t=-251.88e6
            tkill = -251.88e6
            react_enable_bioprod = PB.get_reaction(model, "global", "force_enable_bioprod")
            PB.setvalue!(PB.get_parameter(react_enable_bioprod, "force_times"), [-1e30, tkill, tkill+1.0, 1e30]) 
            PB.setvalue!(PB.get_parameter(react_enable_bioprod, "force_values"), [1.0, 1.0, 0.0, 0.0]) 
        elseif expt == "lowO2"
            # start with low oxygen to test marine sulphur system
            PB.set_variable_attribute!(model, "atm", "O2", :initial_value, 0.1*3.71e19)
        elseif expt == "lowSO4"
            # start with low SO4 to test methane cycling
            PB.set_variable_attribute!!(model, "ocean", "SO4", :initial_value, 100e-3)  # ~100 uM

        elseif expt in ("Sw_T", "Sw_TL", "Sw_2T", "Sw_2Tsx2", "Sw_2Ts")
            # shelf area pH rise mechanism
            shelfarea_force = PB.get_reaction(model, "global", "shelfarea_force")
            reduce_facs = Dict(
                "Sw_T"   => 0.078,              # Take high shelf area case 'Co2HOmL' to low area of 'CoHLOmH'
                "Sw_TL"  => 0.5*0.078,          # Take high shelf area case 'Co2HOmL' to lower than area of 'CoHLOmH'
                "Sw_2T"  => 8.22e10/1.44e12,    # Take high shelf area case 'Co2LOmHW' to low area of 'CoLOmHHW'
                "Sw_2Ts" => 3*8.22e10/1.44e12,  # Take high shelf area case 'Co2LOmHW' to 3*low area of 'CoLOmHHW'
                "Sw_2Tsx2"=>1.5*8.22e10/1.44e12, # Double / half 2Ts
            )
            reduce_fac = reduce_facs[expt]
            PB.setvalue!(PB.get_parameter(shelfarea_force, "force_times"),  [-1e30, -252.05e6, -252.05e6+1.0, 1e30]) 
            PB.setvalue!(PB.get_parameter(shelfarea_force,  "force_values"), [1.0,    1.0,    reduce_fac,   reduce_fac])  

        elseif expt in ("Pp_PE", "Pp_PEe", "Pp_PEes")
            # increase marine prod -> anoxia
            #                  mol P      dur yr  start yr
            Pperts = Dict(
                "Pp_PE"    => (1.3*3e15,   1e5,    -251.95e6 - 2e5),
                "Pp_PEe"   => (1.3*3e15,   2e5,    -251.95e6  -3e5),
                "Pp_PEes"  => (1.3*3e15,   2e5,    -251.95e6  - 3e5)
            )

            (Ptot, duration, tstart) = Pperts[expt]
            tend = tstart + duration
            Ppulse = PB.get_reaction(model, "global", "Ppulse")
            PB.setvalue!(
                PB.get_parameter(Ppulse, "perturb_times"),  
                [-1e30,     tstart,     tstart+1.0,     tend,   tend+1.0, 1e30]
            ) 
            PB.setvalue!(PB.get_parameter(
                Ppulse, "perturb_totals"), 
                [0.0,   0.0,    Ptot/duration,   Ptot/duration,    0.0,             0.0]
            )
            PB.setvalue!(PB.get_parameter(
                Ppulse, "perturb_deltas"), 
                [0.0,   0.0,    0.0,    0.0,    0.0,             0.0]
            )

        elseif expt == "Lk_2"
            # stop land organic carbon burial at 252.0Ma
            force_locbpert = PB.get_reaction(model, "global", "force_locbpert")
            PB.setvalue!(PB.get_parameter(force_locbpert, "force_times"),  [-1e30, -252.0e6, -252.0e6+1.0, 1e30]) 
            PB.setvalue!(PB.get_parameter(force_locbpert, "force_values"), [1.0,    1.0,    0.0,    0.0])  

        elseif expt in ("Cia_1", "Cia_s2")
            # isotopically light carbon injection at 251.95Ma
            CO2pulse_Cia = PB.get_reaction(model, "atm", "CO2pulse_Cia")
            CO2sizes_deltas = Dict(
                "Cia_1"=>(4.86e17, -50.0), 
                "Cia_s2"=>(2/3*2*4.86e17, -25.0)
            )
            (size, delta) = CO2sizes_deltas[expt] # mol C
            duration = 0.5e5 # yr
            PB.setvalue!(PB.get_parameter(CO2pulse_Cia, "perturb_times"),  
                [-1e30, -251.95e6, -251.95e6+1.0, -251.95e6+duration, -251.95e6+duration+1.0, 1e30]) 
            PB.setvalue!(PB.get_parameter(CO2pulse_Cia, "perturb_totals"), 
                [0.0,   0.0,    size/duration,   size/duration,    0.0,             0.0])
            PB.setvalue!(PB.get_parameter(CO2pulse_Cia, "perturb_deltas"),
                delta.*[1.0,   1.0,    1.0,      1.0,            1.0,             1.0])

        elseif expt in ("Cib_1", "Cib_2")
            # isotopically neutral carbon injection at 251.89Ma
            CO2pulse_Cib = PB.get_reaction(model, "atm", "CO2pulse_Cib")
            CO2sizes=Dict("Cib_1"=>2e18, "Cib_2"=>4e18)
            size = CO2sizes[expt]
            ep2d13C = 2.65
            duration = 1e4 # yr
            PB.setvalue!(PB.get_parameter(CO2pulse_Cib, "perturb_times"),  
                [-1e30, -251.89e6, -251.89e6+1.0, -251.89e6+duration, -251.89e6+duration+1.0, 1e30]) 
            PB.setvalue!(PB.get_parameter(CO2pulse_Cib, "perturb_totals"), 
                [0.0,   0.0,    size/duration,   size/duration,    0.0,             0.0])
            PB.setvalue!(PB.get_parameter(CO2pulse_Cib, "perturb_deltas"),
                ep2d13C.*[1.0,   1.0,    1.0,      1.0,            1.0,             1.0])
        
        else
            error("unrecognized expt='$(expt)'")
        end
    end

    return model
end

function plot_PTB3box(
    output;
    pager=PALEOmodel.DefaultPlotPager()
)
    pager(plot(title="Degass Weathering burial", output, 
                ["atm.ccdeg", "land.silw", "land.carbw", "land.oxidw", "fluxOceanBurial.flux_total_Ccarb", 
                 "oceanfloor.shelf_Ccarb_total", "fluxOceanBurial.flux_total_Corg", "land.locb", "atm.Cia", "atm.Cib"], 
                ylabel="flux (mol C yr-1)"))

    return nothing
end
