
using Plots
using Printf

"""
    make_images(output)

Produce single-panel svg figures for ./images

Defaults and titles are for `output` produced by `examples/mitgcm/MITgcm_2deg8_PO4MMcarbSCH4.jl`
"""
function make_images(
    output;
    lon1=198.28, # Pac cell centre for 2.8 deg
    lon2=338.91, # Atl cell centre for 2.8 deg
    toutput = 1e12, # last timestep
    show_sections=true,
    lon_lims = (0.0, 360.0),
    lat_lims = (-90.0, 90.0),
)
    tmodel = last(PALEOmodel.get_array(output, "ocean.tmodel").values)
    println("using tmodel $tmodel (yr)")

    println("longitudes for sections $lon1 (deg), $lon2 (deg)")
       
    println("lon limits $lon_lims, lat limits $lat_lims")
 
    gr(size=(700, 500))
    p = heatmap(
        title=@sprintf("Surface P (mol m-3) at %.2f yr", tmodel), 
        output, "ocean.P_conc", (tmodel=toutput, zt_isel=1, expand_cartesian=true); 
        swap_xy=true, xlabel="lon (deg)", ylabel="lat (deg)", margin=(5, :mm), xlims=lon_lims, ylims=lat_lims,
    )
    if show_sections
        plot!(p, [lon1, lon1],  [-90, 90], linecolor=:red, linewidth=2, label=nothing)
        plot!(p, [lon2, lon2],  [-90, 90], linecolor=:red, linewidth=2, label=nothing)
    end
    display(p)
    savefig(p, "surface_P_2deg8_PO4MMcarbSCH4_2000yr.svg")

    gr(size=(800, 400))
    p = heatmap(
        title=@sprintf("H2S (mol m-3) lon %.2f at %.2f yr", lon1, tmodel),
        output, "ocean.H2S_conc", (tmodel=toutput, lon=lon1, expand_cartesian=true);
        swap_xy=true, xlabel="lat (deg)", ylabel="depth (m)", clims=(0.0, 0.11), margin=(5, :mm),
    )
    display(p)
    savefig(p, "H2S_lon1_2deg8_PO4MMcarbSCH4_2000yr.svg")

    p = heatmap(
        title=@sprintf("H2S (mol m-3) lon %.2f at %.2f yr", lon2, tmodel), 
        output, "ocean.H2S_conc", (tmodel=toutput, lon=lon2, expand_cartesian=true);
        swap_xy=true, xlabel="lat (deg)", ylabel="depth (m)", clims=(0.0, 0.11), margin=(5, :mm),
    ) 
    display(p)
    savefig(p, "H2S_lon2_2deg8_PO4MMcarbSCH4_2000yr.svg")

    return nothing
end