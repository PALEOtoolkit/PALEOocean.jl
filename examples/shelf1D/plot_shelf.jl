
function plot_shelf_phys(
    output;
    pager=PALEOmodel.DefaultPlotPager(),
)
    pager(
        heatmap(title="Temperature (K)", output, "ocean.temp", (column=1,)),
        plot(title="Wind speed",  output, ["oceansurface.wind_speed"], (cell=1,), ylabel="wind speed (m s-1)",),
        heatmap(title="log10(kZ) (m s-1) NB: cell interior faces (z scale off-by-one)", output, "ocean.Kz", (column=1,), map_values=log10),
        plot(title="surface insolation (W m-2)", output, ["oceansurface.surface_insol"], (cell=1,), ylabel="insolation (W m-2)",),
        heatmap(title="opacity (m-1)", output, "ocean.opacity", (column=1,)),
        heatmap(title="insolation (W m-2)", output, "ocean.insol", (column=1,)),
    )
    return nothing
end


function plot_tracers_conc(
    output;
    domain="ocean",
    tracers=[],
    colT=[],
    plot_totals=true,
    pager=PALEOmodel.DefaultPlotPager()
)
    for tr in tracers
        pager(heatmap(title="$domain $tr conc", output, "$(domain).$(tr)_conc", (column=1, )))
        if !isempty(colT)
            pager(plot(title="$domain $tr conc", output, "$(domain).$(tr)_conc", (tmodel=colT, column=1),
                        swap_xy=true, labelattribute=:filter_records))
        end
        if plot_totals
            pager(plot(title="$domain $tr total", output, "$(domain).$(tr)_total"))
        end
    end

    return nothing
end

function plot_tracers_other(
    output;
    domain="ocean",
    tracers=["pHtot", "carbchem/OmegaAR"],
    colT=[],
    pager=PALEOmodel.DefaultPlotPager()
)
    for tr in tracers
        pager(heatmap(title="$domain $tr", output, "$(domain).$(tr)", (column=1, )))
        if !isempty(colT)
            pager(plot(title="$domain $tr", output, "$(domain).$(tr)", (tmodel=colT, column=1),
                        swap_xy=true, labelattribute=:filter_records))
        end
    end

    return nothing
end


function plot_airsea(
    output;
    tracers=["O2"],
    pager=PALEOmodel.DefaultPlotPager,
)
    for tr in tracers
        pager(plot(title="Air-sea $tr", output, ["fluxAtmtoOceansurface.flux_total_$tr"], ylabel="mol yr-1",))
    end
    return nothing
end

function plot_biota(
    output;
    colT=[],
    pager=PALEOmodel.DefaultPlotPager(),
)

    pager(
        plot(title="Total P components", output, ["ocean.P_total", "ocean.DOP_total", "ocean.POP_total", "ocean.phytP_total"], ylabel="P (mol)", areaplot=true),
        plot(title="Total P", output, ["global.total_P"],  ylabel="P (mol)",),
    )
    plot_tracers_conc(output, tracers=["P", "DOP", "POP", "bioprod/phytP"], colT=colT, plot_totals=false, pager=pager)

    pager(
        plot(title="Production Corg", output, ["ocean.Prod_Corg_total"], ylabel="NPP (mol Corg yr-1)",)
    )

    return nothing
end

function plot_biota_S(
    output;
    colT=[],
    pager=PALEOmodel.DefaultPlotPager(),
)
    pager(
        plot(title="Total P components", output, ["ocean.P_total", "ocean.DOP_total", "ocean.POP_total", "ocean.phytP_total", "ocean.phytPS_total"], ylabel="P (mol)", areaplot=true),
        plot(title="Total P", output, ["global.total_P"],  ylabel="P (mol)",),
    )

    plot_tracers_conc(output, tracers=["P", "DOP", "POP", "bioprod/phytP", "bioprodS/phytP"], colT=colT, plot_totals=false, pager=pager)

    pager(
        plot(title="Production Corg", output, ["ocean.Prod_Corg_total", "ocean.ProdS_Corg_total", "ocean.redox_CH4_SO4_total"], ylabel="NPP (mol Corg yr-1)",)
    )

    return nothing
end

function plot_oceanfloor(
    output;
    solutes=["O2", "H2S", "SO4", "CH4"],
    pager=PALEOmodel.DefaultPlotPager,
)
    fluxes = ["fluxOceanfloor.particulateflux_total_Corg"]
    for s in solutes
        push!(fluxes, "fluxOceanfloor.soluteflux_total_"*s)
    end
    pager(
        plot(title="Oceanfloor C, O2, S fluxes", output, fluxes, ylabel="mol yr-1"),
    )

    fluxP = 1/106.0 * PALEOmodel.get_array(output, "fluxOceanfloor.particulateflux_total_Corg")
    pager(
        (
            plot(title="Oceanfloor P fluxes", output, ["fluxOceanfloor.soluteflux_total_P"]);
            plot!(fluxP, ylabel="mol yr-1")
        )

    )

    return nothing
end


function plot_totals_S(
    output;
    pager=PALEOmodel.DefaultPlotPager(),
)
    pager(
        plot(title="Total O2eq", output, ["global.total_O2eq"], ylabel="O2eq (mol)",),
        plot(title="Total S components", output, ["global.total_S", "ocean.H2S_total", "ocean.SO4_total"], ylabel="S (mol)",),
        plot(title="Total S", output, ["global.total_S"], ylabel="S (mol)",),
        plot(title="Total C", output, ["global.total_C"], ylabel="C (mol)",),
        plot(title="Total P", output, ["global.total_P"], ylabel="P (mol)",),
    )

    return nothing
end

