function plot_totals(
    output;
    species=["C", "TAlk", "TAlkerror", "P", "O2", "S"],
    pager=PALEOmodel.DefaultPlotPager(),
    plotargs=NamedTuple(),
)
    # conservation checks
    
    for s in species
        if s == "S"
            pager(plot(title="Total $s", output, ["global.total_$(s)"]; ylabel="$s (mol)",), plotargs...)
        else
            pager(plot(title="Total $s", output, ["ocean.$(s)_total"]; ylabel="$s (mol)",), plotargs...)
        end
    end

    return nothing
end

function plot_airsea(
    output;
    columns=[:-],
    species=["O2", "CO2"],
    pager=PALEOmodel.DefaultPlotPager(),
    plotargs=NamedTuple(),
)
    for sp in species
        pager((plot(title="Air-sea $sp", output, "fluxAtmtoOceansurface.flux_$sp", (cell=columns,); plotargs...);
            plot!(output, "fluxAtmtoOceansurface.flux_total_$sp"; ylabel="$sp flux (mol yr-1)", plotargs...)))
    end

    return nothing
end

function plot_ocean_tracers(
    output;
    columns=[:-],
    colskip=0, # number of blank panels to add (use to format pages)
    tracers=[
        "insol", "DIC_conc", "TAlk_conc", "temp", "pHtot", "O2_conc", "SO4_conc", "H2S_conc", "CH4_conc","P_conc", "Ca_conc",
        "H2S_delta", "SO4_delta", "CH4_delta", "OmegaAR"],
    tcol=[-Inf, 10.0, 100.0, 1000.0, Inf], # model time for column plots
    pager=PALEOmodel.DefaultPlotPager(),
    plotargs=NamedTuple(),
)
    for tr in tracers
        for col in columns
            pager(plot(title="$tr $col", output, "ocean.$tr", ( tmodel=tcol, column=col);
                    swap_xy=true, labelattribute=:filter_records, plotargs...))
        end
        for i in 1:colskip
            pager(:skip)
        end
    end

    return nothing
end

function plot_ocean_tracers_ts(
    output;
    columns=[:-],
    colskip=0, # number of blank panels to add (use to format pages)
    tracers=["O2_conc", "H2S_conc"],
    pager=PALEOmodel.DefaultPlotPager(),
    plotargs=NamedTuple(),
)
    bodgedplotargs = Dict(pairs(plotargs))
    if get(bodgedplotargs, :xflip, false)
        @warn "plot_ocean_tracers_ts: replacing xflip with mult_x_coord=-1.0"
        bodgedplotargs[:mult_x_coord]=-1.0
        delete!(bodgedplotargs, :xflip)
        # fix up xlim if present
        xlim = get(bodgedplotargs, :xlim, nothing)
        if !isnothing(xlim)            
            bodgedplotargs[:xlim] = (-xlim[2], -xlim[1])
        end
    end

    for tr in tracers  
        for col in columns
            pager(
                heatmap(title="Ocean $tr :$col", output, "ocean.$tr", (column=col,); bodgedplotargs...),
            )
        end
        for i in 1:colskip
            pager(:skip)
        end
    end
       
    return nothing
end

