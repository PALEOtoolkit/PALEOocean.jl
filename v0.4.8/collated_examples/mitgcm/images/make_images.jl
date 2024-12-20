
# gr(size=(1200, 900))


using Printf

# produce figures for ./images
function make_images(
    output;
    lonidx1=71, # 2.8 Pac   200 ECCO
    lonidx2=121, # 2.8 Atl  340 ECCO
)
    toutput = 1e12 # last timestep
    tmodel = last(PALEOmodel.get_array(output, "ocean.tmodel").values)
    println("using tmodel $tmodel (yr)")


    # get lon, lat coords from surface P
    P_conc = PALEOmodel.get_array(output, "ocean.P_conc", (k=1, tmodel=toutput))
    lon, lon_lower, lon_upper = P_conc.dims[1].coords
    lat, lat_lower, lat_upper = P_conc.dims[2].coords
    lon1 = lon.values[lonidx1]
    lon2 = lon.values[lonidx2]
    println("longitudes for sections $lonidx1 = $lon1 (deg), $lonidx2 = $lon2 (deg)")    
    # lat_lims = (first(lat_lower.values), last(lat_upper.values))
    lat_lims = (-90.0, 90.0)
    println("lat limits $lat_lims")

 
    gr(size=(700, 500))
    p = heatmap(
        title=@sprintf("Surface P (mol m-3) at %.2f yr", tmodel), output, "ocean.P_conc", (tmodel=toutput, k=1); 
        swap_xy=true, xlabel="lon (deg)", ylabel="lat (deg)", margin=(5, :mm), ylims=lat_lims,
    )
    plot!(p, [lon1, lon1],  [-90, 90], linecolor=:red, linewidth=2, label=nothing)
    plot!(p, [lon2, lon2],  [-90, 90], linecolor=:red, linewidth=2, label=nothing)
    display(p)
    savefig(p, "surface_P_2deg8_PO4MMcarbSCH4_2000yr.svg")

    gr(size=(800, 400))
    p = heatmap(
        title=@sprintf("H2S (mol m-3) lon %.2f at %.2f yr", lon1, tmodel), output, "ocean.H2S_conc", (tmodel=toutput, i=lonidx1);
        swap_xy=true, xlabel="lat (deg)", ylabel="depth (m)", clims=(0.0, 0.11), margin=(5, :mm),
    )
    display(p)
    savefig(p, "H2S_lon1_2deg8_PO4MMcarbSCH4_2000yr.svg")

    p = heatmap(
        title=@sprintf("H2S (mol m-3) lon %.2f at %.2f yr", lon2, tmodel), output, "ocean.H2S_conc", (tmodel=toutput, i=lonidx2);
        swap_xy=true, xlabel="lat (deg)", ylabel="depth (m)", clims=(0.0, 0.11), margin=(5, :mm),
    ) 
    display(p)
    savefig(p, "H2S_lon2_2deg8_PO4MMcarbSCH4_2000yr.svg")

    # heatmap(output, "ocean.$tr", (tmodel=t, i=lonidx2), swap_xy=true), 

    return nothing
end