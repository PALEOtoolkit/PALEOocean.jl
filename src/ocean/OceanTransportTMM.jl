module OceanTransportTMM

import SparseArrays
import MAT   # Matlab file access
import LinearAlgebra
using Printf

import SIMD
import Preferences

# import Infiltrator # Julia debugger

import PALEOboxes as PB
using PALEOboxes.DocStrings
import PALEOocean


"""
    ReactionOceanTransportTMM

GCM ocean transport implementation using transport matrices in format defined by [Khatiwala2007](@cite)
Requires download of Samar Khatiwala's TMM files (for MITgcm, UVic models) as described in
<https://github.com/samarkhatiwala/tmm> where TM files are from <http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/>


NB: `length(operatorID)` must be 2, to define `operatorID[1]` for explicit and `operatorID[2]` for implicit matrices.

The `base_path` parameter sets the top level of the folder structure for the downloaded matrices.

Code based on `TMM/tmm/models/petsc3.4/mitgchem/matlab/make_input_files_for_migchem_dic_biotic_model.m`

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_SETUP)
"""
Base.@kwdef mutable struct ReactionOceanTransportTMM{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParString("base_path", "\$TMMDir\$/MITgcm_2.8deg",
            description="directory containing transport matrices"),

        PB.ParBool("sal_norm", false,
            description="apply salinity normalisation to transport matrix"),

        PB.ParBool("use_annualmean", false,
            description="true to read annual mean matrix"),

        PB.ParInt("num_seasonal", 12, 
            description="number of seasonal matrices"),

        PB.ParInt("Aimp_deltat", 3600*24, units="seconds",
            description="timestep to derive upscaling factor for implicit transport matrix"),

        PB.ParBool("kji_order", true,
            description="true to sort indices into k,j,i order to optimise memory layout"),

        PB.ParInt("pack_chunk_width", 4, allowed_values=[0, 2, 4, 8, 16],
            description="non-zero to enable SIMD packed transport matrix multiply"),

        PB.ParInt("TMfpsize", 64, units="bits", allowed_values=[32, 64],
            description="FP size for transport matrix"),
    )

    index_perm          = nothing  # optional permutation to apply to TM indices to optimise memory layout
    index_perm_inverse  = nothing

    grid_ocean          = nothing
    grid_oceansurface   = nothing
    grid_oceanfloor     = nothing
   

    "Transpose of transport matrices: Units yr^{-1}
     so dc/dt = c * trspt_dtm_tr [yr^-1]
     where c is row vector of tracer concentrations
     
    Multiple matrices are stored in Julia CSC format with a common sparsity pattern"
    trspt_dAexp_tr  = nothing
    trspt_dAimp_tr  = nothing
    Aimp_mult::Int    = -1 # upscaling factor Aimp
    TMeltype        = nothing # set from par_TMfpsize
    pack_datatype   = nothing # set from par_pack_chunk_width

    matrix_times::Vector{Float64} = Vector{Float64}()  # model times for transport matrices
    "time interpolator for transport matrices"
    matrix_tinterp              = nothing

    config_data                 = nothing  # raw config_data from .mat file
    grid                        = nothing
    boxes                       = nothing

    matrix_name::String         = ""
    matrix_path::String         = ""    
    matrix_timestep::Float64    = NaN  
    matrix_start_time::Float64  = NaN
    matrix_delta_time::Float64  = NaN
    matrix_cycle_time::Float64  = NaN
    
end
   
# Provide a custom create_reaction implementation so we can set default operatorID
function PB.create_reaction(::Type{ReactionOceanTransportTMM}, base::PB.ReactionBase)
    rj = ReactionOceanTransportTMM(base=base)
    rj.base.operatorID = [1, 2]
    return rj
end

function PB.set_model_geometry(rj::ReactionOceanTransportTMM, model::PB.Model)

    config_data_path = joinpath(rj.pars.base_path[], "config_data.mat")
    PB.setfrozen!(rj.pars.base_path)

    @info "set_model_geometry $(PB.fullname(rj)) reading config_data from file $(config_data_path)"
    rj.config_data = MAT.matread(config_data_path)

    _read_grids(rj)

    PALEOocean.Ocean.set_model_domains(model, rj.grid_ocean, rj.grid_oceansurface, rj.grid_oceanfloor)    
    
    return nothing
end

function _read_grids(rj::ReactionOceanTransportTMM)
   
    grid_path = joinpath(rj.pars.base_path[], "grid.mat")
    @info "$(PB.fullname(rj)) reading grid from file $(grid_path)"
    rj.grid = MAT.matread(grid_path)

    # NB: we define grid z coords to be -ve
    zmid        = -rj.grid["z"][:,1]
    zupper      = -rj.grid["z"][:,1] .+ rj.grid["dznom"][:,1]./2
    zlower      = -rj.grid["z"][:,1] .- rj.grid["dznom"][:,1]./2
    zedges      = vcat(zupper, zlower[end]) 

    lat         = rj.grid["y"][:,1] # mid-point of each bin ?
    latedges    = Vector(undef, length(lat)+1)
    latedges[1:end-1] = rj.grid["y"][:,1] - rj.grid["dth"][1,:]/2
    latedges[end] = rj.grid["y"][end,1] + rj.grid["dth"][1, end]/2

    lon         = rj.grid["x"][:,1] # mid-point of each bin ?
    lonedges    = Vector(undef, length(lon)+1)
    lonedges[1:end-1] = rj.grid["x"][:,1] - rj.grid["dphi"][:,1]/2
    lonedges[end] = rj.grid["x"][end,1] + rj.grid["dphi"][end,1]/2

    rj.grid_ocean = PB.Grids.CartesianGrid(
        PB.Grids.CartesianLinearGrid,
        ["lon", "lat", "zt"], [length(lon), length(lat), length(zmid)], [lon, lat, zmid], [lonedges, latedges, zedges];
        zdim=3,
        zidxsurface=1,
        ztoheight=1.0
    )     
     
    rj.grid_oceansurface = PB.Grids.CartesianGrid(
        PB.Grids.CartesianLinearGrid,
        ["lon", "lat"], [length(lon), length(lat)], [lon, lat], [lonedges, latedges]
    ) 

    rj.grid_oceanfloor   = deepcopy(rj.grid_oceansurface)

    # read linear cell <-> grid mapping (transport matrix uses linear cell index)
    boxes_path = joinpath(rj.pars.base_path[], rj.config_data["matrixPath"], "Data", "boxes.mat")
    @info "$(PB.fullname(rj)) reading linear cell index <-> 3D grid mapping from file $(boxes_path)"
    rj.boxes = MAT.matread(boxes_path)

    v_i = Int.(rj.boxes["ixBox"][:, 1])
    v_j = Int.(rj.boxes["iyBox"][:, 1])
    v_k = Int.(rj.boxes["izBox"][:, 1])

    if rj.pars.kji_order[]
        @info "  reordering transport matrix indices into k, j, i order"
        (rj.index_perm, rj.index_perm_inverse) = PB.Grids.linear_kji_order(rj.grid_ocean, v_i, v_j, v_k)
        v_i = v_i[rj.index_perm]; v_j = v_j[rj.index_perm]; v_k = v_k[rj.index_perm]
    end
    PB.setfrozen!(rj.pars.kji_order)

    PB.Grids.set_linear_index(rj.grid_ocean, v_i, v_j, v_k)
    @info "  set ocean linear <--> cartesian mapping for $(rj.grid_ocean.ncells) cells"

    lsurf = findall(x->x==rj.grid_ocean.zidxsurface, v_k)
    PB.Grids.set_linear_index(rj.grid_oceansurface, v_i[lsurf], v_j[lsurf])
    PB.Grids.set_linear_index(rj.grid_oceanfloor, v_i[lsurf], v_j[lsurf])
    @info "  set surface, floor linear <--> cartesian mapping for $(rj.grid_oceansurface.ncells) cells"

    # set subdomain mappings
    # println("lsurf ", typeof(lsurf), " =", lsurf)
    # println("v_k ", typeof(v_k), " =", v_k)
    PB.Grids.set_subdomain!(rj.grid_ocean, "oceansurface", PB.Grids.BoundarySubdomain(lsurf), true)
    @info "  set ocean.oceansurface Subdomain size=$(length(lsurf))"
    # find ifloor
    lfloor = similar(lsurf)
    for l in eachindex(lfloor)
        fcart = rj.grid_oceanfloor.cartesian_index[l]
        i,j = fcart[1], fcart[2]
        for k in reverse(1:size(rj.grid_ocean.linear_index, 3))
            if !ismissing(rj.grid_ocean.linear_index[i,j,k])
                lfloor[l] = rj.grid_ocean.linear_index[i,j,k]
                break
            end
        end
    end
    PB.Grids.set_subdomain!(rj.grid_ocean, "oceanfloor", PB.Grids.BoundarySubdomain(lfloor), true)
    @info "  set ocean.oceanfloor Subdomain size=$(length(lfloor))"

    locean = Vector{Union{Missing, Int}}(undef, rj.grid_ocean.ncells)
    fill!(locean, missing)
    locean[lfloor] .= 1:length(lfloor)
    PB.Grids.set_subdomain!(rj.grid_oceanfloor, "ocean", PB.Grids.InteriorSubdomain(locean), true)
    @info "  set oceanfloor.ocean Subdomain size=$(length(locean))"
    
    return nothing
end

function PB.register_methods!(rj::ReactionOceanTransportTMM)

    PB.add_method_setup!(
        rj,
        setup_grid_TMM,
        (PB.VarList_namedtuple(PALEOocean.Ocean.grid_vars_all),)
    )

    PB.add_method_setup!(
        rj,
        setup_transport_TMM, # reads transport matrices
        (),
    )

    return nothing
end

function setup_grid_TMM(
    m::PB.ReactionMethod,
    (grid_vars, ),
    cellrange::PB.AbstractCellRange,
    attribute_name
)
    attribute_name == :setup || return

    rj = m.reaction

    z_coord = rj.grid_ocean.coords[3]
    z_coord_edges = rj.grid_ocean.coords_edges[3]

    da = rj.grid["da"]
    dv = rj.grid["dv"]
    for l in 1:rj.grid_ocean.ncells
        cartidx = rj.grid_ocean.cartesian_index[l]
        kidx = cartidx[3]

        grid_vars.Abox[l]  = da[cartidx]
        grid_vars.zupper[l] = z_coord_edges[kidx]
        grid_vars.zmid[l]   = z_coord[kidx]      
        grid_vars.zlower[l] = z_coord_edges[kidx+1]

        grid_vars.volume[l] = dv[cartidx]
    end
    grid_vars.volume_total[] = sum(grid_vars.volume)

    grid_vars.Asurf .= grid_vars.Abox[rj.grid_ocean.subdomains["oceansurface"].indices]
    grid_vars.Afloor .= grid_vars.Abox[rj.grid_ocean.subdomains["oceanfloor"].indices]
    grid_vars.Afloor_total[] = sum(grid_vars.Afloor)
    grid_vars.zfloor .= grid_vars.zlower[rj.grid_ocean.subdomains["oceanfloor"].indices]

    grid_vars.rho_ref   .= 1027
    grid_vars.pressure  .= -grid_vars.zmid   # pressure(dbar) ~ depth (m)

    return nothing
end     
           


function PB.register_dynamic_methods!(rj::ReactionOceanTransportTMM)

    @info "register_dynamic_methods! $(nameof(typeof(rj))) $(PB.fullname(rj))"
    length(rj.operatorID) == 2 || error("  configuration error: length(operatorID) != 2")

    if rj.pars.TMfpsize[] == 32
        rj.TMeltype = Float32
    elseif rj.pars.TMfpsize[] == 64
        rj.TMeltype = Float64
    else
        error("unknown TMfpsize $(rj.pars.TMfpsize[])")
    end
   
    @info "  setting TMeltype $(rj.TMeltype)"

    vars = [
        PB.VarDepScalar("global.tforce", "yr",  "historical time at which to apply forcings, present = 0 yr"),
    ]

    (transport_conc_vars, transport_sms_vars, _, num_components) =
        PALEOocean.Ocean.find_transport_vars(
            rj.domain, 
            add_transport_input_vars=false,
        )

    if rj.pars.pack_chunk_width[] == 0
        @info "  using unpacked transport"
        rj.TMeltype === Float64 || error("  only Float64 supported")
        # Aexp  operatorID[1]
        PB.add_method_do!(
            rj, 
            do_transport_TMM,
            (   # need to keep to 4 arguments required by prepare_transport
                PB.VarList_namedtuple([PB.VarDep.(PALEOocean.Ocean.grid_vars_ocean); vars]),
                PB.VarList_components(transport_conc_vars),
                PB.VarList_components(transport_sms_vars),
                PB.VarList_nothing(), # not using input_vars
            ),
            preparefn=PALEOocean.Ocean.prepare_transport, # add buffer
            operatorID=[rj.operatorID[1]],
            p=:trspt_dAexp_tr # field name in rj to use
        )
        # Aimp  operatorID[2]
        PB.add_method_do!(
            rj, 
            do_transport_TMM,
            (   # need to keep to 4 arguments required by prepare_transport
                PB.VarList_namedtuple([PB.VarDep.(PALEOocean.Ocean.grid_vars_ocean); vars]),
                PB.VarList_components(transport_conc_vars),
                PB.VarList_components(transport_sms_vars),
                PB.VarList_nothing(), # not using input_vars
            ),
            preparefn=PALEOocean.Ocean.prepare_transport, # add buffer
            operatorID=[rj.operatorID[2]],
            p=:trspt_dAimp_tr # field name in rj to use
        )
    else   
       
        packed_buffer = PALEOocean.Ocean.PackedBuffer(
            num_components,
            rj.grid_ocean.ncells,
            rj.pars.pack_chunk_width[],
            rj.TMeltype
        )
        @info "  packed transport for $num_components concentration components "*
            "with packed_buffer: $(typeof(packed_buffer)) "*
            "pack_eltype $(PALEOocean.Ocean.pack_eltype(packed_buffer)) "*
            "pack_datatype $(PALEOocean.Ocean.pack_datatype(packed_buffer)) "*
            "size packed_conc_array $(size(packed_buffer.packed_conc_array))"

        # dummy Variable to create a dependency to guarantee pack_conc is called before transport
        sequencer_var = PB.VarPropScalar("%reaction%packed_transport_sequencer", "", "dummy Variable to sequence packed transport")
        # pack concentration Variables
        PB.add_method_do!(
            rj, 
            do_transport_TMM_pack_conc,
            (                
                PB.VarList_components(transport_conc_vars),
                PB.VarList_single(sequencer_var),
            ),  
            p=packed_buffer,
        )
        # Aexp  operatorID[1]
        PB.add_method_do!(
            rj, 
            do_transport_TMM_packed,
            (   
                PB.VarList_namedtuple([PB.VarDep.(PALEOocean.Ocean.grid_vars_ocean); vars]),
                PB.VarList_components(transport_sms_vars),
                PB.VarList_single(PB.VarDep(sequencer_var)),
            ),            
            operatorID=[rj.operatorID[1]],
            p=(:trspt_dAexp_tr, packed_buffer),
        )
        # Aimp  operatorID[2]
        PB.add_method_do!(
            rj, 
            do_transport_TMM_packed,
            (   
                PB.VarList_namedtuple([PB.VarDep.(PALEOocean.Ocean.grid_vars_ocean); vars]),
                PB.VarList_components(transport_sms_vars),
                PB.VarList_single(PB.VarDep(sequencer_var)),
            ),
            operatorID=[rj.operatorID[2]],
            p=(:trspt_dAimp_tr, packed_buffer),
        )        

    end

    PB.setfrozen!(rj.pars.TMfpsize, rj.pars.pack_chunk_width)

    return nothing
end


function setup_transport_TMM(
    m::PB.ReactionMethod,
    (),
    cellrange::PB.AbstractCellRange,
    attribute_name,
)
    attribute_name == :setup || return

    rj = m.reaction

    _read_matrix_data(rj)

    return nothing
end

function _read_matrix_data(rj::ReactionOceanTransportTMM)

    # calculate model time for list of matrices
    if rj.pars.use_annualmean[]
        rj.matrix_times = [NaN]
        rj.matrix_tinterp = nothing
    else
        rj.matrix_start_time = 0.5/rj.pars.num_seasonal[]
        rj.matrix_delta_time = 1.0/rj.pars.num_seasonal[]
        rj.matrix_cycle_time = 1.0
        rj.matrix_times = collect(
            range(
                rj.matrix_start_time,
                step=rj.matrix_delta_time, 
                length=rj.pars.num_seasonal[]
            )
        )                
        rj.matrix_tinterp = PB.LinInterp(rj.matrix_times, 1.0)
    end
    PB.setfrozen!(rj.pars.use_annualmean, rj.pars.num_seasonal)

    # transpose, convert datatype, optionally permute
    function permute_indices_transpose(A)
        I, J, V = SparseArrays.findnz(A)
        if rj.pars.kji_order[]
            A = SparseArrays.sparse(rj.index_perm_inverse[J], rj.index_perm_inverse[I], rj.TMeltype.(V)) # swap I, J to transpose
        else
            A = SparseArrays.sparse(J, I, rj.TMeltype.(V))
        end
        return A
    end


    if rj.pars.use_annualmean[]
        Aexp_paths = [joinpath(rj.pars.base_path[], config_data["explicitAnnualMeanMatrixFile"]*".mat")]
        Aexp_name = "Aexpms"
        Aimp_paths = [joinpath(rj.pars.base_path[], config_data["implicitAnnualMeanMatrixFile"]*".mat")]
        Aimp_name = "Aimpms"
        num_matrices = 1
    else
        num_matrices = rj.pars.num_seasonal[]
        Aexp_paths = [
            joinpath(rj.pars.base_path[], rj.config_data["explicitMatrixFileBase"]*@sprintf("_%02i.mat",i))
            for i in 1:num_matrices
        ]
        Aexp_name = "Aexp"
        Aimp_paths = [
            joinpath(rj.pars.base_path[], rj.config_data["implicitMatrixFileBase"]*@sprintf("_%02i.mat",i))
            for i in 1:num_matrices
        ]
        Aimp_name = "Aimp"
    end
            
    # work out upscaling factor for implicit matrix
    deltaT = rj.grid["deltaT"]
    rj.pars.Aimp_deltat[] % deltaT == 0 ||
        error("requested Aimp_deltat $(rj.pars.Aimp_deltat[]) is not a multiple of matrix deltaT = $deltaT")
    rj.Aimp_mult = rj.pars.Aimp_deltat[]/deltaT
    Aimp_deltat_yr = rj.pars.Aimp_deltat[]/PB.Constants.k_secpyr
    PB.setfrozen!(rj.pars.Aimp_deltat)

    # For large (1 deg) matrices, attempt to minimise memory use
    # NB: the underlying hdf5 library is not thread safe so we can't use threads without 
    # adding a lot of complexity / using a lot of memory

    tmp_Aexp_tr = Vector{Any}(nothing, rj.pars.num_seasonal[])
    for i in 1:num_matrices
        Aexp_path = Aexp_paths[i]
        @info "  reading explicit matrix data from $Aexp_path"  
        file = MAT.matopen(Aexp_path)
        matAexp = MAT.read(file, Aexp_name)
        close(file)
        
        tmp_Aexp_tr[i] = permute_indices_transpose(matAexp)
        matAexp = nothing       
        if rj.pars.sal_norm[]
            error("TODO sal_norm")
            # matrix_dt = calcSalNorm(matrix_dt.tocsr(), sal)
        end
        # Aexp is already in differential form, convert s-1 to yr-1
        tmp_Aexp_tr[i] .= PB.Constants.k_secpyr.*tmp_Aexp_tr[i]  
    end
    # convert to CSR format with common sparsity pattern for fast multiply x vector
    rj.trspt_dAexp_tr = PALEOocean.Ocean.create_common_sparsity_tr!(
        tmp_Aexp_tr, 
        do_transpose=false,
        TMeltype=rj.TMeltype
    )
    tmp_Aexp_tr = []

    tmp_Aimp_tr = Vector{Any}(nothing, rj.pars.num_seasonal[])
    for i in 1:num_matrices
        Aimp_path = Aimp_paths[i]
        @info "  reading implicit matrix data from $Aimp_path"  
        file = MAT.matopen(Aimp_path)
        matAimp = MAT.read(file, Aimp_name)
        close(file)

        # convert implicit matrix
        tmp_Aimp_tr[i] = permute_indices_transpose(matAimp)
        matAimp = nothing
        @info "    converting implicit matrix $i deltaT=$deltaT sec x $(rj.Aimp_mult) "*
            "-> $(rj.pars.Aimp_deltat[]) sec ($Aimp_deltat_yr yr)"    
        if rj.pars.sal_norm[]
            error("TODO sal_norm")
            # matrix_dt = calcSalNorm(matrix_dt.tocsr(), sal)
        end
        # Aimp is in 'discrete' form (C_n+1 = Aimp * C_n) with timestep deltaT sec
        # Upscale by factor rj.Aimp_mult to timestep rj.pars.Aimp_deltat[] seconds, Aimp_deltat_yr years,
        # and convert to differential form, units yr-1
        tmp_Aimp_tr[i] = (tmp_Aimp_tr[i]^rj.Aimp_mult - LinearAlgebra.I)./Aimp_deltat_yr
       
    end
    # convert to CSR format with common sparsity pattern for fast multiply x vector
    rj.trspt_dAimp_tr = PALEOocean.Ocean.create_common_sparsity_tr!(
        tmp_Aimp_tr,
        do_transpose=false, 
        TMeltype=rj.TMeltype
    )
    tmp_Aimp_tr = []

    PB.setfrozen!(rj.pars.sal_norm)    

    return nothing
end


function do_transport_TMM(
    m::PB.ReactionMethod,
    (grid_tforce_vars, conc_components, sms_components, _, buffer),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction
    tm_field = m.p

    recswts = PALEOocean.Ocean.get_weights(
        rj.matrix_tinterp,
        grid_tforce_vars.tforce[],
    )
   
    # need to supply a type explicitly as tm member Variable is untyped 
    TrsptCSCType = PALEOocean.Ocean.TrsptCSC{Float64}
    PALEOocean.Ocean.do_transport_tr(
        grid_tforce_vars, conc_components, sms_components, nothing, buffer,
        getfield(rj, tm_field)::TrsptCSCType, recswts,
        cellrange
    )

    return nothing
end


"calculate packed concentration in a separate reaction so all tiles available for do_transport"
function do_transport_TMM_pack_conc(
    m::PB.ReactionMethod, 
    (conc_components, dummy_var),
    cellrange::PB.AbstractCellRange,
    deltat,
)
    packed_buffer = m.p
  
    PALEOocean.Ocean.transport_pack_conc!(
        packed_buffer, conc_components,
        cellrange)

    return nothing
end


function do_transport_TMM_packed(
    m::PB.ReactionMethod,
    (grid_tforce_vars, sms_components, dummy_var),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction
    tm_field, packed_buffer = m.p

    recswts = PALEOocean.Ocean.get_weights(
        rj.matrix_tinterp,
        grid_tforce_vars.tforce[],
    )
   
    # need to supply a type explicitly as tm member Variable is untyped 
    TrsptCSCType = PALEOocean.Ocean.TrsptCSC{eltype(packed_buffer.packed_conc_array)}
    PALEOocean.Ocean.do_transport_tr(
        grid_tforce_vars, sms_components, packed_buffer,    
        getfield(rj, tm_field)::TrsptCSCType, recswts,
        cellrange,
    )

    return nothing
end


"set default TMMDir key if not present in LocalPreferences.toml"
function __init__()    

    defaultTMMdir = normpath(@__DIR__, "../../../TMM")
    for mod in [@__MODULE__, PB]  # forcing Reactions in PALEOboxes may also need TMMDir
        if !Preferences.has_preference(mod, "TMMDir")        
            @info "OceanTransportTMM adding default TMMDir => $defaultTMMdir to [$mod] LocalPreferences.toml (modify for your local setup)"
            Preferences.set_preferences!(mod, "TMMDir" => defaultTMMdir)
        end
    end
    
    return nothing
end

end # module
