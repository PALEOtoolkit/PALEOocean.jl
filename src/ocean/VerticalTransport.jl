"Ocean vertical transport Reactions"
module VerticalTransport

import SparseArrays

import PALEOboxes as PB
using PALEOboxes.DocStrings

# import Infiltrator


"""
    ReactionLightColumn

Calculate light availability `insol` in ocean interior, given surface insolation `surface_insol`.

Includes: (i) a `background_opacity`; (ii) contributions from any Variables representing concentrations with non-initialize_to_zero
`specific_light_extinction` attribute; (iii) any other opacity contributions added to the Target Variable `opacity`.

# Parameters
$(PARS)
"""
Base.@kwdef mutable struct ReactionLightColumn{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("background_opacity", 0.04, units="m-1",
            description="background opacity"),
    )

    "names of Variables with non-zero :specific_light_extinction"
    names_specific_light_extinction::Vector{String} = String[]
    "values read from :specific_light_extinction attributes"
    vec_specific_light_extinction::Vector{Float64} = Float64[] 
end

PB.register_methods!(rj::ReactionLightColumn) = nothing

# find all Variables (representing concentrations) with attribute :specific_light_extinction != 0.0
function _find_opacity_variables(domain::PB.Domain)
    filter_opacity(v) = PB.get_attribute(v, :specific_light_extinction, 0.0) != 0.0
    return PB.get_variables(domain, filter_opacity)
end

function PB.register_dynamic_methods!(rj::ReactionLightColumn)

    vars = [
        PB.VarDepColumn("oceansurface.surface_insol", "W m-2", "surface downwelling radiative flux"),
        PB.VarProp("insol", "W m-2", "interior downwelling radiative flux"),
        PB.VarTarget("opacity", "m-1", "total opacity from all contributions"),

        PB.VarDep("zupper",  "m",    "depth of upper surface of box (m)  0 is surface, -100 is depth of 100 m"),
        PB.VarDep("zlower",  "m",    "depth of lower surface of box (m)  0 is surface, -100 is depth of 100 m"),
        PB.VarDep("zmid",  "m",    "depth of mid point of box (m)"),
    ]

    # find all Variables (representing concentrations) with attribute :specific_light_extinction != 0.0
    # (value will be reread later in setup_light_column to allow configuration updates)   
    rj.names_specific_light_extinction = [v.name for v in _find_opacity_variables(rj.domain)]

    # create ReactionVariables (which will be linked to domvars_conc_opacity)
    vars_conc_opacity = [
        PB.VarDep(name, "", "")
        for name in rj.names_specific_light_extinction
    ]
    io = IOBuffer()
    println(io, "register_dynamic_methods! $(PB.fullname(rj)) "*
        "adding opacity contributions from $(length(vars_conc_opacity)) Variables:")
    for v in vars_conc_opacity
        println(io, "    $(v.localname)")
    end
    @info String(take!(io))

    PB.add_method_setup!(
        rj,
        setup_light_column,
        (PB.VarList_vector(vars_conc_opacity), ),
    )

    PB.add_method_do!(
        rj,
        do_light_column,
        (PB.VarList_namedtuple(vars), PB.VarList_vector(vars_conc_opacity), ),
    )

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end

"read :specific_light_extinction attributes "
function setup_light_column(
    m::PB.ReactionMethod,
    (vars_conc_opacity, ),
    cellrange::PB.AbstractCellRange,
    attribute_name
)
    attribute_name == :setup || return

    rj = m.reaction
    # check no new Variables with non-zero opacity
    current_names = [v.name for v in _find_opacity_variables(rj.domain)]
    new_names = setdiff(current_names, rj.names_specific_light_extinction)
    isempty(new_names) ||
        error("setup_light_column $(PB.fullname(rj)): Variables $new_names have been updated "*
            "from zero to non-zero :specific_light_extinction since model creation "*
            "(fix: set :specific_light_extinction attribute to a dummy but non-zero value in the .yaml config file")

    # read values from Variable :specific_light_extinction attribute
    empty!(rj.vec_specific_light_extinction)
    (vars_conc_opacity, ) = PB.get_variables_tuple(m)
    for v in vars_conc_opacity
        # get specific_light_extinction from the linked VariableDomain
        push!(
            rj.vec_specific_light_extinction,
            PB.get_domvar_attribute(v, :specific_light_extinction)::Float64
        )
    end

    io = IOBuffer()
    println(io, "setup_light_column $(PB.fullname(rj)): Variables with non-zero :specific_light_extinction :")
    for (name, sle) in PB.IteratorUtils.zipstrict(rj.names_specific_light_extinction, rj.vec_specific_light_extinction)
        println(io, "    ", name, sle)
    end
    @info String(take!(io))

    return nothing
end

function do_light_column(
    m::PB.ReactionMethod,
    pars,
    (vars, vars_conc_opacity),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction
   
    @inbounds for (isurf, colindices) in cellrange.columns
        surface_flux = vars.surface_insol[isurf]
        optical_path_length = zero(first(vars.opacity)) # optical path through this column from surface
        for i in colindices
            # Accumulate contributions to opacity in this cell:
            # Background opacity
            vars.opacity[i] += pars.background_opacity[]
            # Contributions from Variables with non-zero :specific_light_extinction :
            for (sle, v_conc) in PB.IteratorUtils.zipstrict(rj.vec_specific_light_extinction, vars_conc_opacity)
                vars.opacity[i] += sle*v_conc[i]
            end

            # Calculate optical path length and light transport:
            # Add contribution from cell upper to cell mid-point
            optical_path_length += vars.opacity[i]*(vars.zupper[i] - vars.zmid[i])
            vars.insol[i] = exp(-optical_path_length)*surface_flux
            # Add contribution from cell mid-point to cell lower (will be used on next iteration)
            optical_path_length += vars.opacity[i]*(vars.zmid[i] - vars.zlower[i])
        end
    end

    return nothing
end


"""
    ReactionExportDirect

Vertical particle sinking represented as instantaneous transport, described by fixed matrices (suitable for small ocean models)

Transports a list of fluxes defined by parameter `fluxlist`, from input Target Variables  with local name
`export_<flux in fluxlist>` to output Contributor Variables `remin_<flux in fluxlist>`.

Matrices are defined by parameter `transportocean` for sinking within water column, and parameter `transportfloor` for flux to ocean floor.

For larger, column-based ocean models use [`ReactionExportDirectColumn`](@ref).
    
# Parameters
$(PARS)
"""
Base.@kwdef mutable struct ReactionExportDirect{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParStringVec("fluxlist", ["P", "N", "Corg"],
            description="names of fluxes to transport"),
        PB.ParDoubleVecVec("transportocean",
            description="matrix describing ocean export"),
        PB.ParDoubleVecVec("transportfloor",
            description="matrix describing oceanfloor export"),

        PB.ParDouble("conserv_errthresh",  1e-5,
            description="error threshold for conservation check (fraction of input)"),
    )

    trspt_ocean_tr::SparseArrays.SparseMatrixCSC{Float64, Int64} = SparseArrays.spzeros(0, 0)
    trspt_floor_tr::SparseArrays.SparseMatrixCSC{Float64, Int64} = SparseArrays.spzeros(0, 0)

    do_transportocean = false
    do_transportfloor = false

    domain_oceanfloor = nothing
end

function PB.register_methods!(rj::ReactionExportDirect, model::PB.Model)

    @info "register_methods! $(PB.fullname(rj)) add ocean fluxlist=$(rj.pars.fluxlist.v)"
    vars_input_target = PB.Fluxes.FluxTarget("export_", rj.pars.fluxlist,
        isotope_data=rj.external_parameters,
        description="input particulate export flux"
    )
    # We need to access vars_input_target from two different methods (ocean and oceanfloor),
    # which isn't possible for a Target variable.
    # So add as a do_nothing Target, then access (twice) as VarDep.
    PB.add_method_do_nothing!(rj, vars_input_target)

    if !isempty(rj.pars.transportocean.v)
        rj.do_transportocean = true
        @info "$(PB.fullname(rj)) add ocean fluxlist=$(rj.pars.fluxlist.v)"
        vars_oceanremin = PB.Fluxes.FluxContrib("remin_", rj.pars.fluxlist,
            isotope_data=rj.external_parameters,
            description="output particulate export flux"
        )
        vars_input = [PB.VarDep(v) for v in vars_input_target]
        PB.add_method_do!(
            rj,
            do_export_direct_ocean,
            (PB.VarList_components(vars_input), PB.VarList_components(vars_oceanremin)),
            # preparefn=prepare_export_direct,
        )
    end

    if !isempty(rj.pars.transportfloor.v)
        rj.do_transportfloor = true
        rj.domain_oceanfloor = PB.get_domain(model, "oceanfloor")
        !isnothing(rj.domain_oceanfloor) ||
            error("$(PB.fullname(rj)) configuration error: no Domain oceanfloor")

        @info "$(PB.fullname(rj)) add oceanfloor fluxlist=$(rj.pars.fluxlist.v)"
        vars_oceanfloor = PB.Fluxes.FluxContrib(
                "fluxOceanfloor.particulateflux_", rj.pars.fluxlist,
                isotope_data=rj.external_parameters,
                description="output particulate export flux"
        )
        # no length check, as we access ocean from oceanfloor Domain
        vars_input_unchecked = []
        for v in vars_input_target
            vv = PB.VarDep(v)
            PB.set_attribute!(vv, :check_length, false; allow_create=true,)
            PB.reset_link_namestr!(vv, "ocean."*v.localname)
            push!(vars_input_unchecked, vv)
        end
        PB.add_method_do!(
            rj,
            do_export_direct_oceanfloor,
            (PB.VarList_components(vars_input_unchecked), PB.VarList_components(vars_oceanfloor)),
            domain=rj.domain_oceanfloor,
        )
    end

    PB.add_method_setup!(
        rj,
        setup_export_direct,
        (),
    )

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end

function setup_export_direct(
    m::PB.ReactionMethod,
    pars,
    _,
    cellrange::PB.AbstractCellRange,
    attribute_name,
)
    rj = m.reaction
    attribute_name == :setup || return

    # calculate matrices and check matrix sizes
    if rj.do_transportocean
        @info "setup_export_direct: $(PB.fullname(rj)) ocean vertical transport from Parameter 'transportocean' = $(pars.transportocean.v)"
           
        trspt_ocean = PB.vecvecpar_matrix(pars.transportocean)
        size(trspt_ocean) == (PB.get_length(rj.domain), PB.get_length(rj.domain)) ||
            error("$(PB.fullname(rj)) invalid Parameter 'transportocean' matrix size != ocean ncells")

        rj.trspt_ocean_tr = SparseArrays.sparse(transpose(trspt_ocean)) # store transpose so it is easy to multiply
    else
        isempty(pars.transportocean.v) ||
            error("$(PB.fullname(rj)) ocean vertical transport not enabled: set Parameter 'transportocean' to non-empty in config file")
    end

    if rj.do_transportfloor
        @info "setup_export_direct: $(PB.fullname(rj)) ocean -> oceanfloor transport from Parameter 'transportfloor' = $(pars.transportfloor.v)"

        trspt_floor = PB.vecvecpar_matrix(pars.transportfloor)
        size(trspt_floor) == (PB.get_length(rj.domain_oceanfloor), PB.get_length(rj.domain)) ||
            error("$(PB.fullname(rj)) invalid Parameter 'transportfloor' matrix size != ocean ncells")

            rj.trspt_floor_tr = SparseArrays.sparse(transpose(trspt_floor)) # store transpose so it is easy to multiply
    else
        isempty(pars.transportfloor.v) || 
            error("$(PB.fullname(rj)) ocean -> oceanfloor transport not enabled: set Parameter 'transportfloor' to non-empty in config file")
    end

    # Check conservation
    for i=1:PB.get_length(rj.domain)
        sumoutput = 0.0
        if rj.do_transportocean
            sumoutput += sum(trspt_ocean[:, i])
        end
        if rj.do_transportfloor
            sumoutput += sum(trspt_floor[:, i])
        end
        abs(sumoutput-1.0) < pars.conserv_errthresh[] ||
            error("Reaction $(PB.fullname(rj)) conservation error exceeds threshold: ocean cell $i output $sumoutput != 1.0")
    end

    return nothing
end


# output[jrange] = input*A_tr[:, jrange]
function _mul_sparse_jrange!(output, jrange, A_tr::SparseArrays.SparseMatrixCSC, input)
    for j in jrange
        for idx in SparseArrays.nzrange(A_tr, j)
            i = A_tr.rowval[idx]
            output[j] += A_tr.nzval[idx]*input[i]
        end
    end
    return nothing
end
# ignore unlinked output variable
_mul_sparse_jrange!(output::Nothing, jrange, A_tr::SparseArrays.SparseMatrixCSC, input) = nothing

function do_export_direct_ocean(
    m::PB.ReactionMethod,
    (components_input, components_oceanremin),
    cellrange::PB.AbstractCellRange,
    deltat
)
    PB.IteratorUtils.check_lengths_equal(
        components_input, components_oceanremin;
        errmsg="do_export_direct_ocean: components length mismatch components_input, components_oceanremin (check :field_data (ScalarData, IsotopeLinear etc) match)"
    )

    for (input, output) in PB.IteratorUtils.zipstrict(components_input, components_oceanremin)
        _mul_sparse_jrange!(output, cellrange.indices, m.reaction.trspt_ocean_tr, input)
    end

    return nothing
end

function do_export_direct_oceanfloor(
    m::PB.ReactionMethod,
    (components_input, components_oceanfloor),
    cellrange::PB.AbstractCellRange,
    deltat
)
    PB.IteratorUtils.check_lengths_equal(
        components_input, components_oceanfloor;
        errmsg="do_export_direct_ocean: components length mismatch components_input, components_oceanfloor (check :field_data (ScalarData, IsotopeLinear etc) match)"
    )

    for (input, output) in PB.IteratorUtils.zipstrict(components_input, components_oceanfloor)
        _mul_sparse_jrange!(output, cellrange.indices, m.reaction.trspt_floor_tr, input)
    end

    return nothing
end


"""
    ReactionExportDirectColumn

Vertical particle sinking represented as instantaneous transport. As [`ReactionExportDirect`](@ref),
but for regular column-based models with functional form of flux vs depth defined by `exportfunction` parameter.

Transports a list of fluxes defined by parameter `fluxlist`, from input Target Variables  with local name
`export_<flux in fluxlist>` to output Contributor Variables `remin_<flux in fluxlist>`.

# Parameters
$(PARS)
"""
Base.@kwdef mutable struct ReactionExportDirectColumn{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParStringVec("fluxlist", String[],
            description="names of fluxes to transport"),

        PB.ParBool("transportfloor", true,
            description="true to provide oceanfloor flux, false to recycle flux into lowest ocean cell"),

        PB.ParString("exportfunction", "SumExp", allowed_values=keys(EXPORT_FUNCTIONS),
            description="functional form for particle flux vs depth"),

        PB.ParDoubleVec("input_frac", [1.0],
            description="fractions of input for each component ([1.0] for Martin, length=number of components for SumExp"),

        PB.ParDoubleVec("sumexp_scale", [500.0], units="m",
            description="length scales for each component of exponential decay of flux with depth"),

        PB.ParDouble("martin_rovera", 0.858,
            description="Martin power law exponent: flux \\propto depth^rovera"),
        PB.ParDouble("martin_depthmin", 100.0, units="m",
            description="Martin power law minimum depth for start of decay with depth"),
    )

    Ncomps::Int64 = -1
end

PB.register_methods!(rj::ReactionExportDirectColumn) = nothing

function PB.register_dynamic_methods!(rj::ReactionExportDirectColumn)

    @info "register_dynamic_methods!: $(PB.fullname(rj)) "*
        "export_function $(rj.pars.exportfunction[]) ocean fluxlist=$(rj.pars.fluxlist.v)"

    vars = [
        PB.VarDep("zupper",  "m",    "depth of upper surface of box (m)  0 is surface, -100 is depth of 100 m"),
        PB.VarDep("zlower",  "m",    "depth of lower surface of box (m)"),
        PB.VarDep("Abox",    "m^2",  "horizontal area of box"),
        # NB: we access Afloor via ocean subdomain which will provide ocean -> oceanfloor indices mapping
        PB.VarDep("oceanfloor.ocean.Afloor",   "m^2",  "horizontal area of seafloor at base of box"),
    ]

    vars_input = PB.Fluxes.FluxTarget("export_", rj.pars.fluxlist,
        isotope_data=rj.external_parameters,
        description="input particulate export flux"
    )

    vars_oceanremin = PB.Fluxes.FluxContrib("remin_", rj.pars.fluxlist,
        isotope_data=rj.external_parameters,
        description="output particulate export flux"
    )

    if rj.pars.transportfloor[]
        @info "    add oceanfloor fluxlist=$(rj.pars.fluxlist.v)"
        vars_oceanfloor = values(PB.Fluxes.FluxContrib(
            "fluxOceanfloor.particulateflux_", rj.pars.fluxlist,
            isotope_data=rj.external_parameters,
            description="output particulate export flux"
        ))
        # no length check here as we access fluxOceanfloor from ocean Domain using indices from Afloor
        vars_oceanfloor_unchecked = PB.set_attribute!.(vars_oceanfloor, :check_length, false; allow_create=true,)
    else
        vars_oceanfloor = []
        vars_oceanfloor_unchecked = vars_oceanfloor
    end

    # export function to use
    export_function, ncomps_function = EXPORT_FUNCTIONS[rj.pars.exportfunction[]]
    PB.setfrozen!(rj.pars.exportfunction)

    PB.add_method_do!(
        rj,
        do_export_direct_column,
        (
            PB.VarList_namedtuple(vars),
            PB.VarList_components(vars_input),
            PB.VarList_components(vars_oceanremin),
            isempty(vars_oceanfloor_unchecked) ? PB.VarList_nothing() : PB.VarList_components(vars_oceanfloor_unchecked),
        ),
        p=(export_function, ncomps_function),
        preparefn=prepare_do_export_direct_column,
    )

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end

function _export_functions_dict()

    # exponential decay of export flux with depth
    function fracfluxout_SumExp(pars, vars, i, ic)
        # fraction of flux into top of box that leaves bottom of box
        frac_fluxout = exp((vars.zlower[i]-vars.zupper[i])/pars.sumexp_scale[ic])
        return frac_fluxout
    end
    # number of components in flux vs depth function
    function ncomps_SumExp(pars, currentNcomps)
        Ncomps = length(pars.sumexp_scale)
        currentNcomps == -1 || currentNcomps == Ncomps ||
            error("ncomps_SumExp $(PB.fullname(rj)): length(sumexp_scale) has changed")
        return Ncomps
    end

    # Martin power-law decay of export flux with depth
    function fracfluxout_Martin(pars, vars, i, ic)
        # fraction of flux into top of box that leaves bottom of box
        zlower_eff = min(vars.zlower[i], -pars.martin_depthmin[])
        zupper_eff = min(vars.zupper[i], -pars.martin_depthmin[])
        frac_fluxout = (zupper_eff/zlower_eff)^pars.martin_rovera[]
        return frac_fluxout
    end
    # number of components in flux vs depth function
    ncomps_Martin(pars, currentNcomps) = 1

    return Dict(
        "SumExp" => (fracfluxout_SumExp, ncomps_SumExp),
        "Martin" => (fracfluxout_Martin, ncomps_Martin),
    )
end

const EXPORT_FUNCTIONS = _export_functions_dict()

# create buffers
function prepare_do_export_direct_column(
    m::PB.ReactionMethod,
    (
        vars,
        components_input,    # Vector of component Arrays
        components_oceanremin,
        components_oceanfloor,
    )
)
    rj = m.reaction
    export_function, ncomps_function = m.p

    rj.Ncomps = ncomps_function(rj.pars, rj.Ncomps)

    (_, rvars_input, _, _) = PB.get_variables_tuple(m)

    # buffer to accumulate per-cell fluxes (one per thread)
    flux_buf = [
        Array{eltype(first(components_input))}(undef, length(components_input), rj.Ncomps)
        for t in 1:Threads.nthreads()
    ]
    # add flux_buf to the end of the Tuple
    return (vars, components_input, components_oceanremin, components_oceanfloor, flux_buf)
end


function do_export_direct_column(
    m::PB.ReactionMethod,
    pars,
    (
        vars,
        components_input,
        components_oceanremin,
        components_oceanfloor,
        flux_buf,
    ),
    cellrange::PB.AbstractCellRange,
    deltat,
)
    rj = m.reaction
    export_function, ncomps_function = m.p

    Ncomps = ncomps_function(pars, rj.Ncomps)    
    PB.check_parameter_sum(pars.input_frac, Ncomps, tol=1e-6) || 
        error("do_export_direct_column: input_frac does not have the correct number of components and sum to 1.0")
    length(components_input) == length(components_oceanremin) || 
        error("do_export_direct_column: components length mismatch components_input, components_oceanremin (check :field_data (ScalarData, IsotopeLinear etc) match)")
    isnothing(components_oceanfloor) || length(components_input) == length(components_oceanfloor) ||
        error("do_export_direct_column: components length mismatch components_input, components_oceanfloor (check :field_data (ScalarData, IsotopeLinear etc) match)")

    flux = flux_buf[Threads.threadid()]

    (data_Afloor, oceanfloor_indices) = vars.Afloor

    for (icol, colindices) in cellrange.columns
        fill!(flux, 0.0)
        for i in colindices
            if ismissing(oceanfloor_indices[i])
                # no oceanfloor under this box
                floor_frac = 0.0
            else
                # calculate areas and partition into oceanfloor flux and sinking flux
                floor_idx = oceanfloor_indices[i]::Int
                floor_frac = data_Afloor[floor_idx]/vars.Abox[i]
            end

            for n in 1:Ncomps
                # fraction of flux into top of box that leaves bottom of box
                fracfluxout = export_function(pars, vars, i, n)
                for j in 1:length(components_input)
                    # flux is flux in to top of cell i
                    flux_remin = (1.0 - fracfluxout)*flux[j,n]
                    components_oceanremin[j][i] += flux_remin

                    flux_out   = fracfluxout*flux[j,n] + components_input[j][i]*pars.input_frac[n]

                    if floor_frac > 0.0
                        floor_flux = floor_frac*flux_out

                        if isnothing(components_oceanfloor)
                            # no explicit ocean floor: assume uniform column and
                            # recycle flux out of bottom box back into bottom box
                            components_oceanremin[j][i] += floor_flux
                        else
                            # flux to ocean floor
                            components_oceanfloor[j][floor_idx] += floor_flux
                        end
                    end

                    flux[j,n] = (1.0 - floor_frac)*flux_out

                end
            end
        end
    end

    return nothing
end

"Simplified version for a regular column with oceanfloor only at base.
Little difference in speed"
function do_export_direct_regular_column(
    m::PB.ReactionMethod,
    pars,
    (
        vars,
        components_input,
        components_oceanremin,
        components_oceanfloor,
        flux_buf,
    ),
    cellrange::PB.AbstractCellRange,
    deltat,
)
    rj = m.reaction
    export_function, ncomps_function = m.p

    Ncomps = ncomps_function(rj, rj.Ncomps)
    PB.check_parameter_sum(rj.pars.input_frac, Ncomps, tol=1e-6) || 
        error("do_export_direct_regular_column: input_frac does not have the correct number of components and sum to 1.0")
    length(components_input) == length(components_oceanremin) || 
        error("do_export_direct_regular_column: components length mismatch components_input, components_oceanremin (check :field_data (ScalarData, IsotopeLinear etc) match)")
    isnothing(components_oceanfloor) || length(components_input) == length(components_oceanfloor) ||
        error("do_export_direct_regular_column: components length mismatch components_input, components_oceanfloor (check :field_data (ScalarData, IsotopeLinear etc) match)")


    flux = flux_buf[Threads.threadid()]

    (data_Afloor, oceanfloor_indices) = vars.Afloor

    for (icol, colindices) in cellrange.columns
        fill!(flux, 0.0)
        lasti = 0
        for i in colindices

            for n in 1:Ncomps
                # fraction of flux into top of box that leaves bottom of box
                fracfluxout = export_function(rj, vars, i, n)
                lastj = 0
                for j in 1:length(components_input)
                    # flux is flux in to top of cell i
                    flux_remin = (1.0 - fracfluxout)*flux[j,n]
                    components_oceanremin[j][i] += flux_remin

                    flux_out   = fracfluxout*flux[j,n] + components_input[j][i]*pars.input_frac[n]

                    flux[j,n] = flux_out
                end
            end
            lasti = i
        end
        if isnothing(components_oceanfloor)
            # no explicit ocean floor: assume uniform column and
            # recycle flux out of bottom box back into bottom box
            for j in 1:length(components_input)
                for n in 1:Ncomps
                    components_oceanremin[j][lasti] += flux[j, n]
                end
            end
        else
            # flux to ocean floor
            floor_idx = oceanfloor_indices[lasti]::Int
            for j in 1:length(components_input)
                for n in 1:Ncomps
                    components_oceanfloor[j][floor_idx] += flux[j, n]
                end
            end
        end

    end

    return nothing
end



"""
    ReactionSinkFloat

Vertical particle advection. 

Applied to all concentration Variables with non-zero attribute `:vertical_movement` (m d-1, +ve upwards),
using naming convention `<name>_conc` to identify `<name>_sms` Deriv Variable to apply to.
An optional variable `<name>_w` may be defined that overrides `:vertical_movement` to define spatially-variable
vertical motion.

# Parameters
$(PARS)
"""
Base.@kwdef mutable struct ReactionSinkFloat{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParBool("transportfloor", true,
            description="true to provide oceanfloor flux, false to recycle flux into lowest ocean cell"),
    )

    var_rootnames::Vector{String} = String[]
    ":vertical_movement attribute"
    var_vertical_movement::Vector{Float64} = Float64[] 
end

PB.register_methods!(rj::ReactionSinkFloat) = nothing

function PB.register_dynamic_methods!(rj::ReactionSinkFloat)

    vars = [
        PB.VarDep("zupper",  "m",    "depth of upper surface of box (m)  0 is surface, -100 is depth of 100 m"),
        PB.VarDep("zlower",  "m",    "depth of lower surface of box (m)"),
        PB.VarDep("Abox",    "m^2",  "horizontal area of box"),
        # NB: we access Afloor via ocean subdomain which will provide ocean -> oceanfloor indices mapping
        PB.VarDep("oceanfloor.ocean.Afloor",   "m^2",  "horizontal area of seafloor at base of box"; attributes=(:check_length=>false,)),
    ]

    rj.var_rootnames = _find_sinkfloat_rootnames(rj.domain)
   
    vars_conc = [PB.VarDep(rn*"_conc", "", "") for rn in rj.var_rootnames]
    vars_sms = [PB.VarContrib(rn*"_sms", "", "") for rn in rj.var_rootnames]
    vars_w = [PB.VarDep("("*rn*"_w)", "", "") for rn in rj.var_rootnames]
    PB.setfrozen!(rj.pars.transportfloor)
    vars_fluxOceanfloor = if rj.pars.transportfloor[]
        [PB.VarContrib("fluxOceanfloor.sinkflux_$(rn)", "", ""; attributes=(:check_length=>false,)) for rn in rj.var_rootnames] 
    else 
        []
    end
       
    PB.add_method_setup!(
        rj,
        setup_sink_float,
        (
            PB.VarList_tuple(vars_conc),
            PB.VarList_tuple(vars_w),
        ),
    )

    PB.add_method_do!(
        rj,
        do_sink_float,
        (
            PB.VarList_namedtuple(vars),
            PB.VarList_tuple(vars_conc),
            PB.VarList_tuple(vars_sms),
            PB.VarList_tuple(vars_w),
            isempty(vars_fluxOceanfloor) ?
                PB.VarList_tuple_nothing(length(vars_conc)) :
                PB.VarList_tuple(vars_fluxOceanfloor),
        ),
    )

    return nothing
end

"Find all variables with attribute :vertical_movement != 0.0, and check name of form <rootname>_conc"
function _find_sinkfloat_rootnames(domain::PB.Domain)

    filter_conc(v) = PB.get_attribute(v, :vertical_movement, 0.0) != 0.0
    domvars_conc_sinkfloat = PB.get_variables(domain, filter_conc)

    rootnames = String[]
    for v in domvars_conc_sinkfloat
        length(v.name) >= 5 && v.name[end-4:end] == "_conc" ||
            error("find_sinkfloat_rootnames: Variable $(PB.fullname(v)) "*
                "has :vertical_movement attribute but is not named _conc")
        push!(rootnames, v.name[1:end-5])       
    end

    return rootnames
end

function setup_sink_float(
    m::PB.ReactionMethod,
    (
        vars_conc,
        vars_w,
    ),
    cellrange::PB.AbstractCellRange,
    attribute_name
)
    attribute_name == :setup || return

    rj = m.reaction
    # rvars_ are VariableReactions corresponding to data arrays vars_
    (rvars_conc, rvars_w) = PB.get_variables_tuple(m)   

    # check for misconfiguration (:vertical_movement changed from 0.0 to non-zero after model creation)
    current_rootnames = _find_sinkfloat_rootnames(rj.domain)
    new_rootnames = setdiff(current_rootnames, rj.var_rootnames)
    isempty(new_rootnames) ||
        error("setup_sink_float $(PB.fullname(rj)): Variables $(new_rootnames.*"_conc") have been updated "*
            "from zero to non-zero :vertical_movement since model creation "*
            "(fix: set :vertical_movement attribute to a dummy but non-zero value in the .yaml config file")
    
    # read :vertical_movement attribute
    empty!(rj.var_vertical_movement)
    io = IOBuffer()
    println(io, "setup_sink_float: $(PB.fullname(rj)):")
    for (var_w, rvar_w, rvar_conc) in PB.IteratorUtils.zipstrict(vars_w, rvars_w, rvars_conc)
        # NB: get :vertical_movement from the Variable we link to, not our local VarDep
        w_val = PB.get_domvar_attribute(rvar_conc, :vertical_movement)
        isa(w_val, Real) ||
            error("setup_sink_float $(PB.fullname(rj)): Variable $(PB.fullname(rvar_conc.linkvar)) :vertical_movement attribute $w_val is not a number")
        push!(rj.var_vertical_movement, w_val) # not used if w_var linked
        if !isnothing(var_w)
            println(io, "    add $(rvar_conc.linkvar.name) using $(rvar_w.linkvar.name) (m d-1)")
        else            
            println(io, "    add $(rvar_conc.linkvar.name) using :vertical_movement $w_val (m d-1)")
        end
    end
    @info String(take!(io))

    return nothing
end

function do_sink_float(
    m::PB.ReactionMethod,
    (
        vars,
        vars_conc,
        vars_sms,
        vars_w,
        vars_fluxOceanfloor,
    ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction

    PB.IteratorUtils.foreach_longtuple_p(
        do_sink_kernel,
        vars_conc, vars_sms, vars_w, rj.var_vertical_movement, vars_fluxOceanfloor,
        (rj, vars, cellrange, deltat)
    )

end

"calculate vertical advection for a single variable"
function do_sink_kernel(
    ac_conc,
    ac_sms,
    ac_w,
    ac_vertical_movement,
    ac_floor,
    (rj, vars, cellrange, deltat)
)

    (data_Afloor, oceanfloor_indices) = vars.Afloor

    @inbounds for (icol, colindices) in cellrange.columns

        for i in colindices
            w = isnothing(ac_w) ? ac_vertical_movement : ac_w[i]

            # convert m d-1 to m yr-1
            w *= PB.Constants.k_daypyr

            # throttle rate for numerical stability
            if deltat > 0.0
                w = sign(w)*min(abs(w), (vars.zupper[i]-vars.zlower[i])/deltat)
            end

            if w > 0 && i > 1
                # upward flux into box above
                flux = ac_conc[i]*w*vars.Abox[i]
                ac_sms[i] -= flux
                ac_sms[i-1] += flux
            elseif w < 0
                # downward flux into box below and oceanfloor
                flux = ac_conc[i]*(-w)*vars.Abox[i]
                ac_sms[i] -= flux

                if ismissing(oceanfloor_indices[i])
                    # no oceanfloor under this box
                    floor_frac = 0.0
                else
                    # calculate areas and partition into oceanfloor flux and sinking flux
                    floor_idx = oceanfloor_indices[i]::Int
                    floor_frac = data_Afloor[floor_idx]/vars.Abox[i]
                    if isnothing(ac_floor)
                        ac_sms[i] += floor_frac*flux # add back oceanfloor flux to current cell
                    else
                        ac_floor[floor_idx] += floor_frac*flux
                    end
                end

                if i == last(colindices)
                    abs(1.0-floor_frac) < 1e-6 || @warn "deepest ocean cell icol $icol not covered by oceanfloor"
                else
                    # transfer flux to cell below
                    ac_sms[i+1] += (1.0-floor_frac)*flux
                end
            end
        end
    end

    return nothing
end

end # module
