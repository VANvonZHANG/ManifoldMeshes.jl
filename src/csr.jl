struct CSRMapping
    offsets::Vector{Int}
    values::Vector{Int}

    function CSRMapping(offsets::Vector{Int}, values::Vector{Int})
        if offsets[1] != 1
            throw(ArgumentError("offsets must start at 1, got $(offsets[1])"))
        end
        if offsets[end] != length(values) + 1
            throw(ArgumentError("last offset mismatch: $(offsets[end]) != $(length(values) + 1)"))
        end
        new(offsets, values)
    end
end

# Uniform constructor (fixed neighbors per entity)
function CSRMapping(n_entities::Int, neighbors_per_entity::Int)
    offsets = [1 + (i - 1) * neighbors_per_entity for i in 1:(n_entities + 1)]
    values = fill(0, n_entities * neighbors_per_entity)
    CSRMapping(offsets, values)
end

# Variable constructor with write pointers (two-pass build)
function CSRMapping(n_entities::Int, neighbor_counts::Vector{Int})
    @assert length(neighbor_counts) == n_entities
    offsets = Vector{Int}(undef, n_entities + 1)
    offsets[1] = 1
    for i in 1:n_entities
        offsets[i + 1] = offsets[i] + neighbor_counts[i]
    end
    values = fill(0, offsets[end] - 1)
    ptrs = copy(offsets[1:(end - 1)])
    return CSRMapping(offsets, values), ptrs
end

Base.length(csr::CSRMapping) = length(csr.offsets) - 1

n_neighbors(csr::CSRMapping, i::Int) = csr.offsets[i + 1] - csr.offsets[i]

function Base.getindex(csr::CSRMapping, i::Int)
    @boundscheck 1 <= i <= length(csr) || throw(BoundsError(csr, i))
    start = csr.offsets[i]
    stop = csr.offsets[i + 1] - 1
    return @view csr.values[start:stop]
end

# Fast fixed-size path (zero allocation, returns NTuple)
function getindex_fixed(csr::CSRMapping, i::Int, ::Val{K}) where {K}
    @boundscheck 1 <= i <= length(csr) || throw(BoundsError(csr, i))
    @boundscheck n_neighbors(csr, i) == K ||
                 throw(BoundsError("CSR row $i has $(n_neighbors(csr, i)) neighbors, expected $K"))
    base = csr.offsets[i] - 1
    ntuple(k -> csr.values[base + k], Val(K))
end

"""
    _permute_uniform_csr(csr::CSRMapping, perm::Vector{Int}, ::Val{K}) where {K}

Permute the values of a uniform CSR mapping (fixed `K` entries per row) according
to `perm`, where `perm[new_id] = old_id`. Offsets are unchanged; only values are
reordered.

Used by HEALPixGrid to apply ring→nested ordering permutation.
"""
function _permute_uniform_csr(csr::CSRMapping, perm::Vector{Int}, ::Val{K}) where {K}
    n = length(csr)
    new_values = Vector{Int}(undef, n * K)
    for i in 1:n
        base_old = csr.offsets[perm[i]] - 1
        base_new = (i - 1) * K
        for j in 1:K
            @inbounds new_values[base_new + j] = csr.values[base_old + j]
        end
    end
    return CSRMapping(csr.offsets[1:(n + 1)], new_values)
end
