function plot_totals(
    output;
    species=["C", "TAlk", "TAlkerror", "P", "O2", "S"],
    pager=PALEOmodel.DefaultPlotPager()
)
    # conservation checks
    if "C" in species
        pager(plot(title="Total C", output, ["global.total_C"], ylabel="C (mol)",))
        pager(plot(title="Total C moldelta", output, ["global.total_C.v_moldelta"], ylabel="C_moldelta (mol * delta)",))
        pager((plot(title="Total C components", output, ["global.total_C", "atm.CO2"]);
               plot!(output, ["ocean.DIC"], (cell=[:s, :h, :d], ), ylabel="C (mol)")))
    end
    if "TAlk" in species
        pager(plot(title="Total TAlk", output, ["ocean.TAlk_total"], ylabel="TAlk (mol)",))
    end
    if "TAlkerror" in species
        pager(plot(title="Ocean TAlk error", output, ["ocean.TAlk_constraint"], (cell=[:s, :h, :d], ), ylabel="TAlk error (mol)",))
    end
    if "P" in species
        pager(plot(title="Total P", output, ["ocean.P_total"], ylabel="P (mol)",))
    end
    if "O2" in species
        pager(plot(title="Total O2", output, ["global.total_O2"], ylabel="O2 (mol)",))
    end
    if "S" in species
        pager(plot(title="Total S", output, ["global.total_S"], ylabel="S (mol)",))
        pager(plot(title="Total S moldelta", output, ["global.total_S.v_moldelta"], ylabel="S_moldelta (mol * delta)",))
    end  

    pager(:newpage)

    return nothing
end

function plot_ocean_tracers(
    output;
    tracers=[
        "DIC_conc", "TAlk_conc", "temp", "pHtot", "O2_conc", "SO4_conc", "H2S_conc", "CH4_conc","P_conc", "Ca_conc",
        "H2S_delta", "SO4_delta", "CH4_delta", "OmegaAR"],
    pager=PALEOmodel.DefaultPlotPager()
)
    for tr in tracers
        pager(plot(title=tr,  output, ["ocean.$tr"], (cell=[:s, :h, :d], ) ))
    end

    return nothing
end


function plot_oaonly_abiotic(
    output;
    pager=PALEOmodel.DefaultPlotPager()
)
    pager((plot(title="d13C", output, ["atm.CO2_delta"]);
          plot!(output, ["ocean.DIC_delta"], (cell=[:s, :h, :d], ), ylabel="per mil",)))
    pager((plot(title="Air-sea CO2", output, ["fluxAtmtoOceansurface.flux_total_CO2"]);
          plot!(output,  "fluxAtmtoOceansurface.flux_CO2", (cell=[:s, :h], ), ylabel="mol yr-1",)))
    pager((plot(title="pCO2", output, ["atm.pCO2atm"], );
          plot!(output, ["ocean.pCO2"], (cell=[:s, :h, :d], ), ylabel="pCO2 (atm)",)))
    if PB.has_variable(output, "atm.pO2PAL")
        pager(plot(title="pO2", output, ["atm.pO2PAL"]))
    end

    return nothing
end

function plot_carb_open(
    output;
    pager=PALEOmodel.DefaultPlotPager()
)
    pager(plot(title="flys (fraction of Ccarb export buried)", output, ["oceanfloor.flys"], (cell=[:s, :h, :d],), ylabel="flys"))
 
    pager((plot(title="Weathering burial", output, ["land.silw", "land.carbw",
                "fluxOceanBurial.flux_total_Ccarb", "oceanfloor.shelf_Ccarb_total", "oceanfloor.deep_Ccarb_total"]);
           plot!(output, "fluxOceanfloor.particulateflux_Ccarb", (cell=:d, ), ylabel="flux (mol C yr-1)", ylim=(0, 3e13))))

    pager((plot(title="Temperature", output, ["global.TEMP"]);
           plot!(output, "ocean.temp", (cell=[:s, :h, :d], ), ylabel="temp (K)",)))
    return nothing
end

function plot_corg_open(
    output;
    pager=PALEOmodel.DefaultPlotPager()
)
    pager(plot(title="P input output", output, 
        ["fluxRtoOcean.flux_P", "fluxOceanBurial.flux_total_P"], ylabel="P flux (mol P yr-1)", ylim=(0, 1e11)))
    pager(plot(title="Corg burial", output,
        ["fluxOceanBurial.flux_total_Corg"], ylabel="Corg burial flux (mol C yr-1)", ylim=(0, 1e13)))
    pager(plot(title="S input output", output,
        ["fluxRtoOcean.flux_SO4", "fluxOceanBurial.flux_total_GYP", "fluxOceanBurial.flux_total_PYR"],
        ylabel="S flux (mol S yr-1)", ylim=(0, 5e12)))

    return nothing
end

