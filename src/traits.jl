# -- TopologyStyle --
abstract type TopologyStyle end
struct IsGrid <: TopologyStyle end
struct IsSemiGrid <: TopologyStyle end
struct IsMesh <: TopologyStyle end

TopologyStyle(::Type{IsGrid}) = IsGrid()
TopologyStyle(::Type{IsSemiGrid}) = IsSemiGrid()
TopologyStyle(::Type{IsMesh}) = IsMesh()
TopologyStyle(::Type{T}) where {T} = IsMesh()
TopologyStyle(m::T) where {T} = TopologyStyle(T)

# -- CellTypeStyle --
abstract type CellTypeStyle end
struct IsUniform{K} <: CellTypeStyle end
struct IsMixed{MAX_K} <: CellTypeStyle end

CellTypeStyle(::Type{T}) where {T} = error("$(T) must implement CellTypeStyle")
CellTypeStyle(m::T) where {T} = CellTypeStyle(T)

# -- PatchStyle --
abstract type PatchStyle end
struct NoPatch <: PatchStyle end
struct MultiPatch{N} <: PatchStyle end

PatchStyle(::Type{T}) where {T} = NoPatch()
PatchStyle(m::T) where {T} = PatchStyle(T)

# -- AbstractLocation --
abstract type AbstractLocation end
struct NodeLoc <: AbstractLocation end
struct CellLoc <: AbstractLocation end
struct EdgeLoc <: AbstractLocation end

# -- MixedCellTopology: stack-allocated variable-length cell topology --
struct MixedCellTopology{MAX_K} <: AbstractVector{Int}
    indices::NTuple{MAX_K, Int}
    len::Int

    function MixedCellTopology{MAX_K}(indices::NTuple{MAX_K, Int}, len::Int) where {MAX_K}
        if len < 0 || len > MAX_K
            throw(ArgumentError("len must be between 0 and MAX_K (=$MAX_K), got $len"))
        end
        return new{MAX_K}(indices, len)
    end
end

function MixedCellTopology(indices::NTuple{MAX_K, Int}, len::Int) where {MAX_K}
    return MixedCellTopology{MAX_K}(indices, len)
end

Base.size(m::MixedCellTopology) = (m.len,)
function Base.getindex(m::MixedCellTopology, i::Int)
    @boundscheck checkbounds(m, i)
    return m.indices[i]
end
Base.IndexStyle(::Type{<:MixedCellTopology}) = IndexLinear()
