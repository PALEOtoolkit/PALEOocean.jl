

import PALEOboxes as PB

import LinearAlgebra
import SparseArrays
import SIMD

# import Infiltrator # julia debugger

########################################
# Grid Variables
#######################################

"""
    Ocean.grid_vars_ocean

Standard ocean variables (properties) provided by a Reaction implementing ocean transport.

NB: contains Property Variables for a Reaction that sets up the grid values.
"""
const grid_vars_ocean = [
    PB.VarPropStateIndep("volume",              "m^3",  "volume of ocean cells"),
    PB.VarPropScalarStateIndep("volume_total",  "m^3",  "total volume of ocean cells"),

    PB.VarPropStateIndep("Abox",                "m^2",  "horizontal area of box"),
 
    # "m^2 (k_nbox,k_nbox) Area(i, j) of horizontal contact between lower surface of box i and box j"
    # Aoverlap::AbstractMatrix = Matrix{Float64}(undef,0,0)

    PB.VarPropStateIndep("zupper",              "m",    "depth of upper surface of box (m)  0 is surface, -100 is depth of 100 m"),
    PB.VarPropStateIndep("zlower",              "m",    "depth of lower surface of box (m)"),
    PB.VarPropStateIndep("zmid",                "m",    "mean depth of box"),


    PB.VarPropStateIndep("rho_ref", "kg m^-3",  "Ocean transport model density conversion factor "*
        "(NB: this an artificial model quantity determined by the offline model for tracer transport, and may be distinct from "*
        "physical ocean density, and from any density anomaly / potl density used by the offline physical model)"),

    PB.VarProp("pressure",           "dbar", "Ocean pressure"),
]

const grid_vars_surfacefloor = [
    # no length check as create and access oceanfloor and oceansurface Variables from ocean Domain
    PB.VarPropStateIndep("oceanfloor.Afloor",   "m^2",  "horizontal area of seafloor at base of box", attributes=(:check_length=>false,)),
    PB.VarPropScalarStateIndep("oceanfloor.Afloor_total",  "m^2",  "total area of seafloor"),
    PB.VarPropStateIndep("oceanfloor.zfloor",   "m",    "depth of ocean floor (m, -ve)", attributes=(:check_length=>false,)),
    PB.VarPropStateIndep("oceansurface.Asurf",  "m^2",  "horizontal area of oceansurface", attributes=(:check_length=>false,)),   
]

const grid_vars_all = vcat(grid_vars_ocean, grid_vars_surfacefloor)

"""
    set_model_domains(model::PB.Model, oceangrid, surfacegrid, floorgrid)

Set model `Domain` sizes and Subdomains for ocean, oceansurface, oceanfloor `Domain`s.
(Helper function for a Reaction implementing ocean transport.)
"""
function set_model_domains(model::PB.Model, oceangrid, surfacegrid, floorgrid)
    @info "Ocean.set_model_domains:"
    oceandom = PB.get_domain(model, "ocean")
    if !isnothing(oceandom)
        oceandom.grid = oceangrid
        @info "  ocean Domain size=$(oceangrid.ncells) grid=$oceangrid"
    else
        error("No ocean domain found")
    end

    for domname in ["fluxAtmtoOceansurface", "oceansurface"]
        dom = PB.get_domain(model, domname)
        if !isnothing(dom)
            dom.grid = surfacegrid
            @info "  $(domname) Domain size=$(surfacegrid.ncells) grid=$surfacegrid"          
        end
    end

    for domname in ["fluxOceanfloor", "fluxOceanBurial", "oceanfloor"]
        dom = PB.get_domain(model, domname) 
        if !isnothing(dom)
            dom.grid = floorgrid
            @info "  $(domname) Domain size=$(floorgrid.ncells) grid=$floorgrid"         
        end
    end 

    return nothing
end


"""
    find_transport_vars(domain::PB.AbstractDomain; transport_input_components=false) 
        -> (conc_vars, sms_vars, input_components, num_components)

Find all variables in `domain` with attribute :advect == true, and then use naming convention `<rootname>_conc` to
identify `<rootname>_sms` Variables to add transport flux to.

If `transport_input_components = true`, also define `input_components` as `<varname>_transport_input`
to calculate advective transport input into each cell (slow)
"""
function find_transport_vars(
    domain::PB.AbstractDomain;
    add_transport_input_vars::Bool=false)

    # find all Variables in domain with :advect attribute
    filter_conc(v) = PB.get_attribute(v, :advect, false)
    conc_domvars_advect = PB.get_variables(domain, filter_conc)

    conc_vars = PB.VarDepT[]
    sms_vars = PB.VarContribT[]
    if add_transport_input_vars
        input_vars = PB.VarPropT[]
    else
        input_vars = nothing
    end

    # count number of components
    num_components = 0
   
    for v in conc_domvars_advect        
        v.name[end-4:end] == "_conc" ||
            error("find_transport_vars: Variable $(PB.fullname(v)) has :advect attribute == true but is not named _conc")
        # TODO define eg a :totalname attribute instead of a naming convention
        rootname = v.name[1:end-5]
      
        # add local Variables to link to _conc and _sms (matching on name)
        # default :field_data to link to any ScalarData, IsotopeLinear etc
        push!(sms_vars, PB.VarContrib(rootname*"_sms", "mol yr-1", ""))
        push!(conc_vars, PB.VarDep(rootname*"_conc", "mol m-3", ""))
    
        # optionally, calculate transport input into each cell
        # VarProp so need to set :field_data
        if add_transport_input_vars
            push!(input_vars, PB.VarProp(rootname*"_transport_input", "mol yr-1", "", attributes=(:field_data=>PB.get_attribute(v, :field_data),)))
        end

        # count components
        num_components += PB.num_components(v)
    end
   
    @info "find_transport_vars Domain $(domain.name) length(transport_vars)=$(length(sms_vars)) num_components=$num_components"

    return (conc_vars, sms_vars, input_vars, num_components)
end





