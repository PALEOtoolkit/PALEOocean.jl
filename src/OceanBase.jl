

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

##########################################
# Transport matrices
##########################################

"""
    add_loop!(A::AbstractMatrix, vm::AbstractVector, L::Real, loopindices)

Add circulation with flux `L` to transport matrix `A`, around a closed loop `loopindices`.

Transport matrix `A` (s^{-1}) represents tracer transport,
     
    dc/dt = A * c
     
where `c` is Vector of tracer concentrations

`loopindices` is a list of cell indices representing a closed loop, eg [2, 3, 4, 2]

Units for tracers and fluxes are:

- Ocean / volume based:
    - `c`:  mol m-3 tracer concentration
    - `vm`: m^3. volume per cell
    - `L`:  m^3 s-1 (volume flux, cf 1 Sverdrup = 1e6 m^3 s-1)
- Atmosphere / mass based
    - `c`:  kg / total kg, tracer mass mixing ratio 
    - `vm`: kg, total mass per call
    - `L`:  kg s-1 (mass flux)
"""
function add_loop!(A::AbstractMatrix, vm::AbstractVector, L::Real, loopindices)
    first(loopindices) == last(loopindices) || 
        error("add_loop!: loopindices $loopindices is not a closed loop")

    for (i1, i2) in PB.IteratorUtils.zipstrict(loopindices[begin:end-1], loopindices[begin+1:end])
        A[i2, i1] += L/vm[i2]
        A[i1, i1] -= L/vm[i1]
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
        push!(sms_vars, PB.VarContrib(rootname*"_sms", "", ""))
        push!(conc_vars, PB.VarDep(rootname*"_conc", "", ""))
    
        # optionally, calculate transport input into each cell
        # VarProp so need to set :field_data
        if add_transport_input_vars
            push!(input_vars, PB.VarProp(rootname*"_transport_input", "", "", attributes=(:field_data=>PB.get_attribute(v, :field_data),)))
        end

        # count components
        num_components += PB.num_components(v)
    end
   
    @info "find_transport_vars Domain $(domain.name) length(transport_vars)=$(length(sms_vars)) num_components=$num_components"

    return (conc_vars, sms_vars, input_vars, num_components)
end


"""
    prepare_transport(m::ReactionMethod, (grid_vars, conc_components, sms_components, input_components)) -> 
        (grid_vars, conc_components, sms_components, input_components, buffer)

Add an additional `buffer` for do_transport. `conc_components`, `sms_components`, optional `input_components` are
Vectors of data Arrays (as created by VarList_components from lists generated by [`find_transport_vars`](@ref)).
"""
function prepare_transport(m::PB.ReactionMethod, (grid_vars, conc_components, sms_components, input_components))

    if isempty(conc_components)
        buffer = nothing
    else
        buffer = [similar(first(conc_components)) for t in 1:Threads.nthreads()]        
    end

   
    return (grid_vars, conc_components, sms_components, input_components, buffer)
end

"""
    do_transport(grid_vars, conc_components, sms_components, input_components, buffer, dtm::AbstractMatrix, cr::AbstractCellRange)

Calculate transport rates. 
# Parameters
- `grid_vars, conc_components, sms_components, input_components`:  Vectors of data Arrays (created by VarList_components)
- `buffer`: buffer arrays created by [`prepare_transport`](@ref)
- `dtm::AbstractMatrix`: transport matrix (units yr-1)
"""
function do_transport(
    grid_vars, conc_components, sms_components, input_components::Nothing, buffer, 
    dtm::AbstractMatrix, 
    cr::PB.AbstractCellRange
)

    # transport only
    errmsg = "do_transport: components length mismatch conc_components, sms_components (check :field_data (ScalarData, IsotopeLinear etc) match)"
    for (conc_component, sms_component) in PB.IteratorUtils.zipstrict(conc_components, sms_components; errmsg=errmsg)            
        @inbounds for j in 1:size(dtm, 2)
            for i in cr.indices            
                sms_component[i] += dtm[i, j]*conc_component[j]*grid_vars.volume[i]
            end
        end
    end

    return nothing
end

function do_transport(
    grid_vars, conc_components, sms_components, input_components, buffer, 
    dtm::AbstractMatrix, 
    cr::PB.AbstractCellRange
)
 
    # transport and input
    errmsg = "do_transport: components length mismatch conc_components, sms_components, input_components (check :field_data (ScalarData, IsotopeLinear etc) match)"
    for (conc_component, sms_component, input_component) in PB.IteratorUtils.zipstrict(conc_components, sms_components, input_components; errmsg=errmsg)
        @inbounds for i in cr.indices
            input_component[i] = 0.0
        end

        @inbounds for j in 1:size(dtm, 2)
            for i in cr.indices
                t = dtm[i, j]*conc_component[j]*grid_vars.volume[i]          
                sms_component[i] += t
                # we want (transport without the diagonal 'output' terms)*conc
                if i != j
                    input_component[i] += t
                end
            end
        end
    end

    return nothing
end

function do_transport(
    grid_vars, conc_components, sms_components,
    dtm::LinearAlgebra.Tridiagonal, cr::PB.AbstractCellRange
)
    # transport only
    errmsg = "do_transport: components length mismatch conc_components, sms_components (check :field_data (ScalarData, IsotopeLinear etc) match)"
    for (conc_component, sms_component) in PB.IteratorUtils.zipstrict(conc_components, sms_components; errmsg=errmsg)
        @inbounds for i in cr.indices
            sms_component[i] += dtm.d[i]*conc_component[i]*grid_vars.volume[i]
            if i > 1
                sms_component[i] += dtm.dl[i-1]*conc_component[i-1]*grid_vars.volume[i] 
            end
            if i < length(dtm.d)
                sms_component[i] += dtm.du[i]*conc_component[i+1]*grid_vars.volume[i]
            end
        end
    end
  
    return nothing
end

"Memory-bandwidth optimised version of do_transport using transpose of dtm.
This is an optimisation specifically tied to the Compressed Sparse Column storage layout (Julia SparseMatrixCSC)"
function do_transport_tr(
    grid_vars, conc_components, sms_components, input_components, buffer, 
    dtm_tr::SparseArrays.SparseMatrixCSC, cr::PB.AbstractCellRange
)
    length(conc_components) == length(sms_components) || 
        error("do_transport_tr: components length mismatch conc_components, sms_components (check :field_data (ScalarData, IsotopeLinear etc) match)")

    isempty(sms_components) && return nothing # no transport needed (early return to avoid use of unallocated buffer etc)

    do_transport_tr(
        grid_vars, conc_components, sms_components, input_components, buffer, 
        dtm_tr.colptr, dtm_tr.rowval, (dtm_tr.nzval, dtm_tr.nzval), (1.0, 0.0),
        cr
    )

    return nothing
end


"""
    TrsptCSC

Store multiple SparseArrays.SparseMatrixCSC with a common sparsity pattern.

SparseMatrixCSC format, except with a Vector of nzval, ie:
- `colptr`: Column j is in colptr[j]:(colptr[j+1]-1)
- `rowval`: Row indices of stored values
- `nzval`: Vector of Vector of non-zero values
"""
struct TrsptCSC{T}
    colptr::Vector{Int64}
    rowval::Vector{Int64}
    nzval::Vector{Vector{T}}
end

# steady-state - just use a single matrix
get_weights(tinterp::Nothing, tforce) = ((1.0, 1), (0.0, 0))
# linearly interpolate
get_weights(tinterp::PB.LinInterp, tforce) = PB.interp(tinterp, PB.value_ad(tforce))

function do_transport_tr(
    grid_vars, conc_components, sms_components, input_components, buffer, 
    dtm_trs::TrsptCSC, ((wt1, rec1), (wt2, rec2)),
    cr::PB.AbstractCellRange
)
    length(conc_components) == length(sms_components) || 
        error("do_transport_tr: components length mismatch conc_components, sms_components (check :field_data (ScalarData, IsotopeLinear etc) match)")

    isempty(sms_components) && return nothing # no transport needed (early return to avoid use of unallocated buffer etc)

    do_transport_tr(
        grid_vars, conc_components, sms_components, input_components, buffer, 
        dtm_trs.colptr, dtm_trs.rowval, (dtm_trs.nzval[rec1], dtm_trs.nzval[rec2]), (wt1, wt2),
        cr
    )

    return nothing
end


"Memory-bandwidth optimised version of do_transport using transpose of dtm, linearly interpolating two matrices
This is an optimisation specifically tied to the Compressed Sparse Column storage layout (Julia SparseMatrixCSC)"
function do_transport_tr(
    grid_vars, conc_components, sms_components, input_components::Nothing, buffer,
    colptr, rowval, (nzval1, nzval2), (wt1, wt2), 
    cr::PB.AbstractCellRange
)
    PB.IteratorUtils.check_lengths_equal(
        conc_components, sms_components;
        errmsg="do_transport_tr: components length mismatch conc_components, sms_components (check :field_data (ScalarData, IsotopeLinear etc) match)"
    )

    tbuffer = buffer[Threads.threadid()]

    # Memory bandwidth optimisation: apply each column of dtm_tr to the whole set of vectors,
    # then move on to the next column, etc
    @inbounds for jcol in cr.indices
        # precalculate this column of dtm_tr
        bidx = 0
        for idx in colptr[jcol]:colptr[jcol+1]-1  # column jcol of dtm_tr
            irow = rowval[idx]
            bidx += 1
            tbuffer[bidx] = (nzval1[idx]*wt1 + nzval2[idx]*wt2)*grid_vars.volume[jcol]
        end
        # apply to whole set of vectors
        for (conc_component, sms_component) in PB.IteratorUtils.zipstrict(conc_components, sms_components)                                
            bidx = 0
            for idx in colptr[jcol]:colptr[jcol+1]-1  # column jcol of dtm_tr
                irow = rowval[idx]
                bidx += 1 
                # sms_components[cr.indices]' = (conc_components' * dtm_tr[:, cr.indices]).*volume[cr.indices]'
                sms_component[jcol] += conc_component[irow]*tbuffer[bidx]
            end               
        end
    end
  
    return nothing
end


function do_transport_tr(
    grid_vars, conc_components, sms_components, input_components, buffer, 
    colptr, rowval, (nzval1, nzval2), (wt1, wt2),
    cr::PB.AbstractCellRange
)
    PB.IteratorUtils.check_lengths_equal(
        conc_components, sms_components, input_components;
        errmsg="do_transport_tr: components length mismatch conc_components, sms_components, input_components (check :field_data (ScalarData, IsotopeLinear etc) match)"
    )

    tbuffer = buffer[Threads.threadid()]

    # track transport input into each cell
    
    # Memory bandwidth optimisation: apply each column of dtm_tr to the whole set of vectors,
    # then move on to the next column, etc
    @inbounds for jcol in cr.indices
        # precalculate this column of dtm_tr
        bidx = 0
        for idx in colptr[jcol]:colptr[jcol+1]-1  # column jcol of dtm_tr
            irow = rowval[idx]
            bidx += 1
            tbuffer[bidx] = (nzval1[idx]*wt1 + nzval2[idx]*wt2)*grid_vars.volume[jcol]
        end
        # apply to whole set of vectors
        for (conc_component, sms_component, input_component) in PB.IteratorUtils.zipstrict(conc_components, sms_components, input_components)
            input_component[jcol] = 0.0
            bidx = 0
            for idx in colptr[jcol]:colptr[jcol+1]-1  # column jcol of dtm_tr
                irow = rowval[idx]
                bidx += 1
                # sms_component[cr.indices]' = (conc_component' * dtm_tr[:, cr.indices]).*volume[cr.indices]'
                s = conc_component[irow]*tbuffer[bidx]
                sms_component[jcol] += s
                # we want (transport without the diagonal 'output' terms)*conc
                if jcol != irow
                    input_component[jcol] += s
                end
            end             
        end
    end
 
    return nothing
end

"""
    PackedBuffer

Buffer Arrays for optimized SIMD transport using elements of type SIMD.Vec{pack_chunk_width, pack_eltype}.
Creates `packed_conc_array` for packed and padded Variable concentrations, and an additional workspace `buffer`.
Uses additional type parameters to allow `do_transport_tr` specialization on number of Variable components
etc as well as SIMD chunk width.
"""
struct PackedBuffer{pack_eltype, pack_datatype, pack_chunk_width, num_components, num_chunks, packed_conc_pad_width}
    packed_conc_array::Array{pack_eltype, 2}   # packed and padded Variable concentrations
    packed_conc_vec::Vector{pack_eltype}       # reshaped packed_conc_array
    buffer::Vector{Vector{pack_datatype}}      # workspace for transport

    function PackedBuffer(
        num_components,  # number of Variable components to transport
        ncells,           # number of ocean cells
        pack_chunk_width, # SIMD width (eg 4)
        pack_eltype,      # SIMD element type (eg Float64)
    )

        pack_datatype = SIMD.Vec{pack_chunk_width, pack_eltype}

        num_chunks = Int(ceil(num_components/pack_chunk_width)) # number of SIMD-width chunks needed for num_components
        packed_conc_pad_width = num_chunks*pack_chunk_width   # padded number of components (rounding up to SIMD width)
        packed_conc_array = zeros(pack_eltype, packed_conc_pad_width, ncells)
        packed_conc_vec = vec(packed_conc_array)
           
        buffer = [Vector{pack_datatype}(undef, num_chunks) for t in 1:Threads.nthreads()]
        
        return new{pack_eltype, pack_datatype, pack_chunk_width, num_components, num_chunks, packed_conc_pad_width}(
            packed_conc_array,
            packed_conc_vec,
            buffer,
        )
        
    end
end

pack_eltype(pb::PackedBuffer) = eltype(pb.packed_conc_array)
pack_datatype(pb::PackedBuffer) = eltype(eltype(pb.buffer))

"pack concentration Variables into a contiguous array for optimized SIMD transport"
function transport_pack_conc!(packed_buffer::PackedBuffer, conc_components, cr::PB.AbstractCellRange)
    packed_conc_array = packed_buffer.packed_conc_array
    @inbounds for i in cr.indices
        for p in 1:length(conc_components)
            packed_conc_array[p, i] = conc_components[p][i]
        end
    end
    return nothing
end

"Optimized SIMD transport matrix x contiguous packed conc data array"
function do_transport_tr(
    grid_vars, sms_components, packed_buffer::PackedBuffer,    
    dtm_trs::TrsptCSC, ((wt1, rec1), (wt2, rec2)),
    cr::PB.AbstractCellRange
)
    isempty(sms_components) && return nothing # no transport needed (early return to avoid use of unallocated buffer etc)

    kernel_transport_tr_packed(
        packed_buffer,
        grid_vars.volume,
        sms_components,
        dtm_trs.colptr, dtm_trs.rowval, (dtm_trs.nzval[rec1], dtm_trs.nzval[rec2]), (wt1, wt2),
        cr
    )
 
    return nothing
end


"Optimized SIMD transport matrix x contiguous packed conc data array kernel. 
 Uses type parameters to specialize on number of components and SIMD chunk width"
function kernel_transport_tr_packed(
    packed_buffer::PackedBuffer{
        pack_eltype, pack_datatype, pack_chunk_width, num_components, num_chunks, packed_conc_pad_width
    },    
    volume, 
    sms_components,
    colptr, rowval, (nzval1, nzval2), (wt1, wt2), 
    cr
) where {pack_eltype, pack_datatype, pack_chunk_width, num_components, num_chunks, packed_conc_pad_width}

    num_components == length(sms_components) || 
        error("kernel_transport_tr_packed: components length mismatch num_components, sms_components (check :field_data (ScalarData, IsotopeLinear etc) match)")

    packed_conc_vec = packed_buffer.packed_conc_vec
    pack_sms_buffer = packed_buffer.buffer[Threads.threadid()]

    # Memory bandwidth optimisation: apply each column of dtm_tr to the whole set of vectors,
    # then move on to the next column, etc
    @inbounds for jcol in cr.indices     
        for cidx in 1:num_chunks   
            pack_sms_buffer[cidx] = zero(pack_datatype)
        end      
        for idx in colptr[jcol]:colptr[jcol+1]-1  # column jcol of dtm_tr        
            # tm = pack_DataType((nzval1[idx]*wt1 + nzval2[idx]*wt2)*volume[jcol])
            tm = (nzval1[idx]*wt1 + nzval2[idx]*wt2)*volume[jcol]
            irow = rowval[idx]
            conc_vec_idx = packed_conc_pad_width*(irow-1) + 1
            for cidx in 1:num_chunks        
                pack_sms_buffer[cidx] += SIMD.vload(pack_datatype, packed_conc_vec, conc_vec_idx)*tm
                conc_vec_idx += pack_chunk_width                
            end
        end
        cidx = 1
        for pidx in 1:pack_chunk_width:num_components
            for p in pidx:min(pidx+pack_chunk_width-1, num_components)             
                sms_components[p][jcol] += pack_sms_buffer[cidx][p-pidx+1]
            end
            cidx += 1
        end
    end

    return nothing
end


"""
    create_common_sparsity_tr!(a_matrices; do_transpose, TMeltype=Float64) -> TrsptTr

Reduce a collection of matrices to common sparsity pattern, optionally transposing.

NB: `a_matrices` is used as workspace and contents deleted.
"""
function create_common_sparsity_tr!(a_matrices; do_transpose, TMeltype=Float64)    
    @info "  create_common_sparsity_tr! adding $(length(a_matrices)) matrices eltype $(eltype(a_matrices[1])) -> $TMeltype"

    different_sparsity = false
    a_1 = a_matrices[1]
    for i =  2:length(a_matrices)
        a = a_matrices[i]
        if (a.colptr != a_1.colptr) ||  (a.rowval != a_1.rowval)
            different_sparsity = true
            @info "   sparsity different"
            break
        end
    end

    if do_transpose
        @info "  calculating transpose"    
        for i in eachindex(a_matrices)
            a_matrices[i] = SparseArrays.sparse(transpose(a_matrices[i])) # NB: transpose will create a transpose wrapper, not transpose to CSC format       
        end
    end

    tr_nzval = Vector{Vector{TMeltype}}()
    if different_sparsity
        # add transpose of matrices to find common sparsity pattern
        a_n = size(first(a_matrices), 1)
        a_tr_sum = SparseArrays.spzeros(a_n, a_n)
   
        for i in eachindex(a_matrices)      
            a_tr_sum += a_matrices[i]
        end
    
        # sparsity pattern of sum is the common sparsity pattern
        tr_colptr = a_tr_sum.colptr
        tr_rowval  = a_tr_sum.rowval

        # get I,J,V for a_sum
        Isum, Jsum, Vzeros = SparseArrays.findnz(a_tr_sum)
        fill!(Vzeros, 0.0)

        @info "  create_common_sparsity_tr converting to common sparsity pattern nnz=$(length(Vzeros))"
        for i in eachindex(a_matrices)
            I, J, V = SparseArrays.findnz(a_matrices[i])
            a_matrices[i] = SparseArrays.spzeros(eltype(a_matrices[i]), 0, 0) # free memory
            # reconstruct matrix with zeros for each (non-zero) element in a_sum
            # this will store zeros if necessary to match the sparsity pattern of a_sum
            tmp_a_tr = SparseArrays.sparse(append!(I, Isum), append!(J, Jsum), append!(TMeltype.(V), TMeltype.(Vzeros)), 
                            size(a_tr_sum,1 ), size(a_tr_sum,2))
            push!(tr_nzval, tmp_a_tr.nzval)
            tr_colptr == tmp_a_tr.colptr || error("colptr sparsity pattern error")
            tr_rowval == tmp_a_tr.rowval || error("rowval sparsity pattern error")
            length(tmp_a_tr.nzval) == length(Vzeros) || error("nzval sparsity pattern error")
        end
        @info "  create_common_sparsity_tr done"
    else
        # already  common sparsity pattern
        @info "  create_common_sparsity_tr already common sparsity pattern nnz=$(length(a_matrices[1].nzval))"
        tr_colptr = a_matrices[1].colptr
        tr_rowval  = a_matrices[1].rowval
        
        for i in eachindex(a_matrices)
            push!(tr_nzval, TMeltype.(a_matrices[i].nzval))
            a_matrices[i] = SparseArrays.spzeros(eltype(a_matrices[i]), 0, 0) # free memory
        end
    end

    return TrsptCSC(tr_colptr, tr_rowval, tr_nzval)
end
