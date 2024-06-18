function plot_totals(
    output;
    species=["T", "P", "O2", "S"],
    pager=PALEOmodel.DefaultPlotPager(),
    plotargs=NamedTuple(),
)
    # conservation checks
    
    for sp in species
        if PB.has_variable(output, "global.total_$sp")
            pager(plot(title="Total $sp (global)", output, ["global.total_$sp"]; ylabel="$sp (mol)",), plotargs...)
        else # ocean only
            pager(plot(title="Total $sp (ocean)", output, ["ocean.$(sp)_total"]; ylabel="$sp (mol)",), plotargs...)        
        end
    end

    return nothing
end

function plot_airsea(
    output;
    columns=[:hlat, :gyre, :upw],
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
    columns=[:hlat, :gyre, :upw],
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
    columns=[:gyre, :upw],
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

function plot_ocean_floor(
    output;
    columns=[:hlat, :gyre, :upw],
    colskip=0,
    tracers=["flys"],
    tcol=1e12,
    pager=PALEOmodel.DefaultPlotPager(),
    plotargs=NamedTuple(),
)
    for tr in tracers
        for col in columns
            p = plot(title="$tr $col"; plotargs...)
            for t in tcol
                t_a = PALEOmodel.get_array(output, "oceanfloor.$tr", tmodel=t, column=col,)
                # add zdepth dimension 
                zdepth = PB.NamedDimension("zdepth", PALEOmodel.get_array(output, "oceanfloor.zfloor", tmodel=t, column=col,).values)
                t_a = PB.FieldArray(tr, t_a.values, (zdepth, ), t_a.attributes)
                plot!(title="$tr $col", t_a; swap_xy=true, labelattribute=:filter_records, plotargs...)
            end
            pager(p)
        end
        for i in 1:colskip
            pager(:skip)
        end
    end

    return nothing
end

function plot_ocean_burial_CorgP(
    output;
    pager=PALEOmodel.DefaultPlotPager(),
    plotargs=NamedTuple(),
)
    
    pager(
        plot(title="P input output", output,  ["fluxRtoOcean.flux_P", "fluxOceanBurial.flux_total_P"];
            ylabel="P flux (mol P yr-1)", plotargs...),
        plot(title="Production burial", output,  ["ocean.Prod_Corg_total", "fluxOceanBurial.flux_total_Corg"],
            ylabel="C flux (mol C yr-1)", plotargs...),
        (
            plot(title="O2 Corg burial", output,  ["fluxOceanBurial.flux_total_Corg"],
                ylabel="O2, C flux (mol yr-1)", ylim=(-5e12, 5e12), plotargs...);
            plot!(-1*PALEOmodel.get_array(output, "atm.O2_restoring"));
            plot!(-1*PALEOmodel.get_array(output,  "fluxAtmtoOceansurface.flux_total_O2"));
        )
   )
    
    return nothing
end

function plot_ocean_burial(
    output;
    ocean_only=true,  # extra restoring flux contributions
    pager=PALEOmodel.DefaultPlotPager(),
    plotargs=NamedTuple(),
)
   
    pager(
        plot(title="P input output", output,  ["fluxRtoOcean.flux_P", "fluxOceanBurial.flux_total_P"];
             ylabel="P flux (mol P yr-1)", plotargs...),
        plot(title="Production burial", output,  ["ocean.Prod_Corg_total", "ocean.Prod_Ccarb_total", "fluxOceanBurial.flux_total_Corg", "oceanfloor.deep_Ccarb_total"];
             ylabel="C flux (mol C yr-1)", plotargs...),
        plot(title="S input output", output,  ["fluxRtoOcean.flux_SO4", "fluxOceanBurial.flux_total_GYP", "fluxOceanBurial.flux_total_PYR"];
            ylabel="S flux (mol S yr-1)", plotargs...),
        plot(title="TAlk input", output,  ["fluxRtoOcean.flux_TAlk"];
            ylabel="flux (mol yr-1)", plotargs...),
    )
    
    if ocean_only
        pager(
            plot(title="Corg Pyr burial O2 flux", output,  ["fluxOceanBurial.flux_total_Corg", "fluxOceanBurial.flux_total_PYR",
                "atm.O2_restoring", "fluxAtmtoOceansurface.flux_total_O2"];
                ylabel="flux (mol yr-1)", plotargs...),
            plot(title="C input output", output,  ["fluxOceanBurial.flux_total_Corg", "fluxOceanBurial.flux_total_Ccarb",
                "oceanfloor.shelf_Ccarb_total", "oceanfloor.deep_Ccarb_total", "oceanfloor.sfw_total", "atm.CO2_restoring"];
                ylabel="flux (mol yr-1)", plotargs...),
        )
    else
        pager(
            plot(title="Corg Pyr burial O2 flux", output,  ["fluxOceanBurial.flux_total_Corg", "fluxOceanBurial.flux_total_PYR", "fluxAtmtoOceansurface.flux_total_O2"];
                ylabel="flux (mol yr-1)", ylim=(-1e13, 1e13), plotargs...),
            plot(title="C input output", output,  ["fluxOceanBurial.flux_total_Corg", "fluxOceanBurial.flux_total_Ccarb",
                "oceanfloor.shelf_Ccarb_total", "oceanfloor.deep_Ccarb_total", "oceanfloor.sfw_total"];
                ylabel="flux (mol yr-1)", plotargs...),
        )
    end

    return nothing
end



