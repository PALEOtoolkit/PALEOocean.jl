
function plot_forcings(
    output; 
    pager=PALEOmodel.DefaultPlotPager(),
    lonidx=71, # 2.8 Pac   200 ECCO
)

    for t in [0.0, 0.5]
        pager(
            heatmap(output, "oceansurface.wind_speed", (tmodel=t,), swap_xy=true),
            heatmap(output, "oceansurface.open_area_fraction", (tmodel=t,), swap_xy=true),
        )     
    end
    for t in [0.0, 0.5]
        pager(heatmap(output, "oceansurface.surface_insol", (tmodel=t,), swap_xy=true))
    end

    pager(
        heatmap(output, "ocean.insol", (tmodel=0.0, i=lonidx), swap_xy=true),
        heatmap(output, "ocean.insol", (tmodel=0.0, i=lonidx), swap_xy=true, ylim=(-500.0, 0)), 
    )

    for t in [0.0, 0.5]
        pager(
            heatmap(output, "ocean.temp", (tmodel=t, k=1), swap_xy=true),
        )
    end
     
    return nothing
end

function plot_abiotic_O2(
    output;
    toutputs=[0.0, 1.0, 10.0, 100.0, 1000.0],
    pager=PALEOmodel.DefaultPlotPager(),
    lonidx1=71, # 2.8 Pac   200 ECCO
    lonidx2=121, # 2.8 Atl  340 ECCO
)
    for t in toutputs
        pager(
            heatmap(output, "ocean.O2_conc", (tmodel=t, i=lonidx1), swap_xy=true),
            heatmap(output, "ocean.O2_conc", (tmodel=t, i=lonidx2), swap_xy=true), 
        )
    end
    
    pager(
        plot(title="O2 air-sea flux",  output,  ["fluxAtmtoOceansurface.flux_total_O2"], ylabel="flux (mol yr-1)",),
        plot(title="O2 totals", output, ["global.total_O2", "atm.O2", "ocean.O2_total"], ylabel="total (mol)",),
        plot(title="O2 total", output, ["global.total_O2"], ylabel="total (mol)",),
    )
    return nothing
end

function plot_PO4MMbase(
    output;
    toutputs=[0.0, 1.0, 10.0, 100.0, 1000.0],
    tbioprod=[1999.5, 2000.0],
    pager=PALEOmodel.DefaultPlotPager(),
)

    for t in toutputs
        pager(
            heatmap(output, "ocean.P_conc", (tmodel=t, k=1), swap_xy=true,),
        )
    end

    if PB.has_variable(output, "ocean.DOP_conc")
        pager(
            heatmap(output, "ocean.DOP_conc", (tmodel=toutputs[end], k=1), swap_xy=true,),
            plot(title="P total", output, ["global.total_P", "ocean.P_total", "ocean.DOP_total"], ylabel="total (mol)",),
        )
    end

    if PB.has_variable(output, "global.total_P")
        pager(
            plot(title="P total", output, ["global.total_P"], ylabel="total (mol)",),
        )
    elseif PB.has_variable(output, "ocean.P_total")
        pager(
            plot(title="P total", output, ["ocean.P_total"], ylabel="total (mol)",),
        )
    end

    for t in tbioprod
        pager(
            heatmap(output, "ocean.bioprod/Prod_Corg", (tmodel=t, k=1), swap_xy=true,),
            heatmap(output, "ocean.bioprod/Prod_Corg", (tmodel=t, k=2), swap_xy=true,),
        )
    end

    pager(plot(title="Marine production", output, ["ocean.Prod_Ccarb_total", "ocean.Prod_Corg_total"], 
        ylabel="Production (mol C yr-1)",))

    return nothing
end


function plot_carbSCH4(
    output;
    pager=PALEOmodel.DefaultPlotPager(),
)

    pager(
        plot(title="O2 totals", output, ["atm.O2", "ocean.O2_total"], ylabel="total (mol)",),
        plot(title="O2eq total", output, ["global.total_O2eq"], ylabel="total (mol)",),
        plot(title="P total", output, ["global.total_P", "ocean.P_total"], ylabel="total (mol)",),
        plot(title="C total", output, ["global.total_C"], ylabel="total (mol)",),
        plot(title="C total moldelta", output, ["global.total_C.v_moldelta"], ylabel="mol S * per mil",),
        plot(title="S totals", output, ["global.total_S", "ocean.SO4_total", "ocean.H2S_total"], ylabel="total (mol)",),
        plot(title="S total", output, ["global.total_S"], ylabel="total (mol)",),
        plot(title="S total moldelta", output, ["global.total_S.v_moldelta"], ylabel="mol S * per mil",),
        plot(title="TAlk total", output, ["ocean.TAlk_total"], ylabel="total (mol)",),
    
        plot(title="Atmospheric d13CO2", output, ["atm.CO2_delta"], ylabel="d13C (per mil)",),
        plot(title="Atmospheric pCO2", output, ["atm.pCO2atm"], ylabel="pCO2 (atm)",),
        plot(title="Atmospheric pO2", output, ["atm.pO2PAL"], ylabel="pO2 (PAL)",),
    )

    return nothing
end

function plot_tracers(
    output;
    tracers=["SO4_conc", "H2S_conc", "CH4_conc", "SO4_delta", "H2S_delta", "CH4_delta", "TAlk_conc", "DIC_conc", "DIC_delta"],
    toutputs=[1e12],
    pager=PALEOmodel.DefaultPlotPager(),
    lonidx1=71, # 2.8 Pac   200 ECCO
    lonidx2=121, # 2.8 Atl  340 ECCO
)

    for tr in tracers
        for t in toutputs
            pager(
                heatmap(output, "ocean.$tr", (tmodel=t, i=lonidx1), swap_xy=true),
                heatmap(output, "ocean.$tr", (tmodel=t, i=lonidx2), swap_xy=true), 
                heatmap(output, "ocean.$tr", (tmodel=t, k=1), swap_xy=true),
                :skip,
            )
        end
    end

    return nothing
end
