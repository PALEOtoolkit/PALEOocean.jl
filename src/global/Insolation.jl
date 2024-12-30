module Insolation


import PALEOboxes as PB
using PALEOboxes.DocStrings


"""
    ReactionForceInsolationModernEarth
 
Calculate time and latitude dependent daily mean modern Earth surface solar insolation.

Daily mean photosynthetically-active surface insolation 

    `insolation` = TOA flux * (1 - `albedo`) * `parfrac`

See [`insolMITgcmDIC`](@ref) for details.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionForceInsolationModernEarth{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("albedo", 0.6,
            description="mean planetary albedo"),
        PB.ParDouble("parfrac", 1.0,
            description="fraction of radiation that is photosynthetically active"),
        PB.ParDoubleVec("latitude", Float64[], units="degrees N",
            description="if non-empty, override grid latitude and set explicitly for each surface cell"),
    )
end

function PB.register_methods!(rj::ReactionForceInsolationModernEarth)

    vars = [
        PB.VarDepScalar("global.tforce", "yr",  "historical time at which to apply forcings, present = 0 yr"),
        PB.VarProp("insolation", "W m-2", "daily mean surface insolation"),
    ]
    if isempty(rj.pars.latitude.v)
        push!(vars, 
            PB.VarDepScalar("lat",  "",  "coordinate variable"; 
                attributes=(:field_data=>PB.ArrayScalarData, :data_dims=>("lat",),))
        )
    end

    PB.add_method_do!(
        rj,
        do_force_insolation, 
        (PB.VarList_namedtuple(vars),),
        p=rj.domain.grid,  # add as context so fully typed
    )
end


function PB.check_configuration(rj::ReactionForceInsolationModernEarth, model::PB.Model)
    configok = true
    if !isempty(rj.pars.latitude) && !isnothing(rj.domain.grid)
        if rj.domain.grid.ncells != length(rj.pars.latitude)
            @warn "check_configuration $(PB.fullname(rj))  length(latitude) parameter $(length(rj.pars.latitude)) != grid.ncells $(grid.ncells)"
            configok = false
        end
    end
    return configok
end

function do_force_insolation(m::PB.ReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)
    
    grid = m.p
    
    tforce = PB.value_ad(vars.tforce[])

    for i in cellrange.indices
        if isempty(pars.latitude)
            lat = vars.lat[PB.Grids.get_lat_idx(grid, i)]
        else
            lat = pars.latitude[i]
        end
        vars.insolation[i] = insolMITgcmDIC(tforce, lat, albedo=pars.albedo[], parfrac=pars.parfrac[])
    end

    return nothing
end

"""
    insolMITgcmDIC(Timeyr,latdeg; albedo=0.6, solar=1360.0, parfrac=1.0) -> sfac

MITgcm DIC package insol function directly translated from fortran.
Similar to [Brock1981](@cite).

NB: there are three normalization constants here: `solar`, `albedo`, `parfrac` to define
top-of-atmosphere flux (from astronomical formulae) -> a crude approx to ground level flux (taking into account clouds etc) -> photosynthetic PAR flux

    C !DESCRIPTION:
    C find light as function of date and latitude
    C based on paltridge and parson

# Arguments:
- `Timeyr`: yr, model time, NB: year assumed to start in winter
- `latdeg`: deg, latitudes
- `albedo`: planetary albedo (ie correct for top-of-atmosphere to ground-level, clouds etc)
- `solar`: W m-2 solar constant
- `parfrac`: photosynthetically active fraction

# Returns:
- `sfac`: daily average photosynthetically active solar radiation just below surface
"""
function insolMITgcmDIC(Timeyr,latdeg; albedo=0.6, solar=1360.0, parfrac=1.0)
  
    # C find day (****NOTE for year starting in winter*****)
    dayfrac= Timeyr - floor(Timeyr) # fraction of year
    yday = 2*π*dayfrac                    # convert to radians
    delta = (0.006918                
        -(0.399912 *cos(yday))           # cosine zenith angle
        +(0.070257 *sin(yday))           # (paltridge+platt)
        -(0.006758 *cos(2*yday))             
        +(0.000907 *sin(2*yday))
        -(0.002697 *cos(3*yday))
        +(0.001480 *sin(3*yday)) )

    # C latitude in radians
    lat=latdeg*2*π/360

    sun1 = -sin(delta)/cos(delta) * sin(lat)/cos(lat)
    sun1 = max(sun1,-0.999)
    sun1 = min(sun1, 0.999 )
    dayhrs = abs(acos(sun1))
    cosz = ( sin(delta)*sin(lat)+             # average zenith angle
         (cos(delta)*cos(lat).*sin(dayhrs)./dayhrs) )
    cosz = max(cosz, 5e-3)
    frac = dayhrs/π                           # fraction of daylight in day
    # C daily average photosynthetically active solar radiation just below surface
    fluxi = solar*(1-albedo)*cosz*frac*parfrac

    # C convert to sfac
    sfac = max(1e-5,fluxi)

    return sfac
end



end