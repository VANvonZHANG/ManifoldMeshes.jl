# -- TopologyStyle --

"""
    TopologyStyle

Abstract type for mesh topology classification.
Dispatch on `TopologyStyle(mesh)` to get `IsGrid()`, `IsSemiGrid()`, or `IsMesh()`.
"""
abstract type TopologyStyle end

"""
    IsGrid

Trait for fully regular structured grids (e.g., LatLonGrid).
Neighbors follow simple Cartesian (i, j) indexing.
"""
struct IsGrid <: TopologyStyle end

"""
    IsSemiGrid <: TopologyStyle

Trait for structured grids with specialized indexing (e.g., reduced Gaussian,
HEALPix) where simple Cartesian neighbor lookup does not apply.
"""
struct IsSemiGrid <: TopologyStyle end
"""
    IsMesh

Trait for fully unstructured meshes with arbitrary connectivity.
"""
struct IsMesh <: TopologyStyle end

TopologyStyle(::Type{IsGrid}) = IsGrid()
TopologyStyle(::Type{IsSemiGrid}) = IsSemiGrid()
TopologyStyle(::Type{IsMesh}) = IsMesh()
TopologyStyle(::Type{T}) where {T} = IsMesh()
TopologyStyle(m::T) where {T} = TopologyStyle(T)

# -- CellTypeStyle --

"""
    CellTypeStyle

Trait abstract type for cell topology classification.
Concrete subtypes: `IsUniform{K}`, `IsMixed{MAX_K}`.
"""
abstract type CellTypeStyle end

"""
    IsUniform{K} <: CellTypeStyle

Trait indicating all cells have exactly K nodes (e.g., `IsUniform{4}()` for
quadrilateral meshes, `IsUniform{3}()` for triangular meshes).
"""
struct IsUniform{K} <: CellTypeStyle end

"""
    IsMixed{MAX_K} <: CellTypeStyle

Trait indicating cells have varying node counts, capped at MAX_K.
"""
struct IsMixed{MAX_K} <: CellTypeStyle end

CellTypeStyle(::Type{T}) where {T} = error("$(T) must implement CellTypeStyle")
CellTypeStyle(m::T) where {T} = CellTypeStyle(T)

# -- PatchStyle --

"""
    PatchStyle

Trait abstract type for patch/domain structure classification.
Concrete subtypes: `NoPatch`, `MultiPatch{N}`.
"""
abstract type PatchStyle end

"""
    NoPatch <: PatchStyle

Trait indicating a single contiguous domain with no patch subdivision.
"""
struct NoPatch <: PatchStyle end

"""
    MultiPatch{N} <: PatchStyle

Trait indicating N independent structured patches (e.g., 6 faces of a cubed
sphere). Use with `cell_face` and `cell_local_2d` for patch-aware queries.
"""
struct MultiPatch{N} <: PatchStyle end

PatchStyle(::Type{T}) where {T} = NoPatch()
PatchStyle(m::T) where {T} = PatchStyle(T)

# -- ProjectionStyle --

"""
    ProjectionStyle

Abstract type for cubed-sphere projection classification.
Concrete subtypes: `Gnomomic`, `Equiangular`.

A projection changes cell geometry on the manifold, so it is represented as a
type parameter on `CubedSphereGrid` (not a runtime field). By contrast, a
variant that only changes numbering — e.g. HEALPix `:ring` vs `:nested`, where
the cells are identical — stays a runtime `Symbol` field.
"""
abstract type ProjectionStyle end

"""
    Gnomomic <: ProjectionStyle

Gnomonic (tangent-plane) cubed-sphere projection.
"""
struct Gnomomic <: ProjectionStyle end

"""
    Equiangular <: ProjectionStyle

Equiangular cubed-sphere projection.
"""
struct Equiangular <: ProjectionStyle end

function ProjectionStyle(::Type{T}) where {T}
    error("$(T) does not have a ProjectionStyle (only CubedSphereGrid uses one)")
end
ProjectionStyle(m::T) where {T} = ProjectionStyle(T)

# -- AbstractLocation --

"""
    AbstractLocation

Abstract type for staggered-grid location tags.
Used to indicate where data lives relative to the mesh: nodes, cell centers, or edge midpoints.
"""
abstract type AbstractLocation end

"""
    NodeLoc

Location tag: data lives at mesh nodes (vertices).
"""
struct NodeLoc <: AbstractLocation end

"""
    CellLoc

Location tag: data lives at cell centers.
"""
struct CellLoc <: AbstractLocation end

"""
    EdgeLoc

Location tag: data lives at edge midpoints.
"""
struct EdgeLoc <: AbstractLocation end

# -- MixedCellTopology: stack-allocated variable-length cell topology --

"""
    MixedCellTopology{MAX_K} <: AbstractVector{Int}

Stack-allocated variable-length cell topology. Stores up to `MAX_K` integer
indices, with active length `len` (0 ≤ len ≤ MAX_K).

# Examples
```jldoctest
julia> m = MixedCellTopology((1, 2, 3, 0, 0, 0), 3)
3-element MixedCellTopology{6}:
 1
 2
 3
```
"""
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

# Tuple comparison: MixedCellTopology is a value type; compare elementwise
# against plain tuples so `cell_nodes(g, c) == (1, 4, 3)` reads naturally.
function Base.:(==)(m::MixedCellTopology, t::Tuple{Vararg{Int}})
    m.len == length(t) && all(i -> m.indices[i] == t[i], 1:(m.len))
end
Base.:(==)(t::Tuple{Vararg{Int}}, m::MixedCellTopology) = m == t
