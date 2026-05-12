"""
    AbstractManifoldMesh{M <: AbstractManifold}

Abstract supertype for all mesh types on a manifold `M`.

Every concrete subtype must implement the full mesh interface:
grid properties (`manifold`, `num_cells`, `num_nodes`, `num_edges`),
geometry queries (`node_coordinates`, `cell_volume`, `cell_centroid`),
topology queries (`cell_nodes`, `cell_cells`, `node_cells`, `cell_edges`),
edge information (`edge_length`, `edge_midpoint`, `edge_outward_normal`),
and boundary markers (`boundary_nodes`, `boundary_edges`).

Cell, node, and edge IDs are 1-indexed integers.
"""
abstract type AbstractManifoldMesh{M <: AbstractManifold} end

# -- Grid Properties --

"""
    manifold(g) -> AbstractManifold

Returns the underlying manifold on which the mesh is defined.
"""
function manifold(g::AbstractManifoldMesh)
    error("$(typeof(g)) must implement `manifold`")
end

"""
    num_cells(g) -> Int

Total number of cells in the mesh.
"""
function num_cells(g::AbstractManifoldMesh)
    error("$(typeof(g)) must implement `num_cells`")
end

"""
    num_nodes(g) -> Int

Total number of nodes (vertices) in the mesh.
"""
function num_nodes(g::AbstractManifoldMesh)
    error("$(typeof(g)) must implement `num_nodes`")
end

"""
    num_edges(g) -> Int

Total number of edges in the mesh.
"""
function num_edges(g::AbstractManifoldMesh)
    error("$(typeof(g)) must implement `num_edges`")
end

# -- Geometry Queries --

"""
    node_coordinates(g, node_id) -> SVector{D, Float64}

Cartesian coordinates of node `node_id` in the embedding space.
"""
function node_coordinates(g::AbstractManifoldMesh, node_id::Int)
    error("$(typeof(g)) must implement `node_coordinates`")
end

"""
    cell_volume(g, cell_id) -> Float64

Volume (area on 2D manifolds) of cell `cell_id`, pre-computed at construction.
"""
function cell_volume(g::AbstractManifoldMesh, cell_id::Int)
    error("$(typeof(g)) must implement `cell_volume`")
end

"""
    cell_centroid(g, cell_id) -> SVector{D, Float64}

Centroid of cell `cell_id` on the manifold, pre-computed at construction.
"""
function cell_centroid(g::AbstractManifoldMesh, cell_id::Int)
    error("$(typeof(g)) must implement `cell_centroid`")
end

# -- Topology Queries --

"""
    cell_nodes(g, cell_id) -> NTuple{K, Int}

Node IDs at the corners of cell `cell_id`. Tuple length depends on cell type.
For quadrilateral cells, returns (SW, SE, NE, NW).
"""
function cell_nodes(g::AbstractManifoldMesh, cell_id::Int)
    error("$(typeof(g)) must implement `cell_nodes`")
end

"""
    cell_cells(g, cell_id) -> NTuple{K, Int}

Cell IDs of cells sharing a face with `cell_id`. Uses `0` as sentinel for missing neighbors (e.g., at domain boundaries).
For quadrilateral cells, returns (south, north, west, east).
"""
function cell_cells(g::AbstractManifoldMesh, cell_id::Int)
    error("$(typeof(g)) must implement `cell_cells`")
end

"""
    node_cells(g, node_id) -> Vector{Int}

Cell IDs of all cells adjacent to `node_id`.

Note: Allocates a `Vector{Int}` per call. Not recommended for tight FVM loops.
Future unstructured meshes should consider a CSR layout for allocation-free access.
"""
function node_cells(g::AbstractManifoldMesh, node_id::Int)
    error("$(typeof(g)) must implement `node_cells`")
end

"""
    cell_edges(g, cell_id) -> NTuple{K, Int}

Edge IDs of the edges bounding cell `cell_id`.
For quadrilateral cells, returns (south, north, west, east).
"""
function cell_edges(g::AbstractManifoldMesh, cell_id::Int)
    error("$(typeof(g)) must implement `cell_edges`")
end

# -- Edge Information --

"""
    edge_length(g, edge_id) -> Float64

Geodesic arc length of edge `edge_id` on the manifold.
"""
function edge_length(g::AbstractManifoldMesh, edge_id::Int)
    error("$(typeof(g)) must implement `edge_length`")
end

"""
    edge_midpoint(g, edge_id) -> SVector{D, Float64}

Geodesic midpoint of edge `edge_id` on the manifold.
For degenerate (zero-length) edges, returns one of the endpoints.
"""
function edge_midpoint(g::AbstractManifoldMesh, edge_id::Int)
    error("$(typeof(g)) must implement `edge_midpoint`")
end

"""
    edge_outward_normal(g, edge_id, cell_id) -> NamedTuple{(:base_point, :normal)}

Outward-pointing unit normal at the midpoint of `edge_id`, relative to `cell_id`.

Returns `(base_point, normal)` where:
- `base_point` is the midpoint on the manifold (tangent space origin)
- `normal` is the unit normal vector in the tangent space at `base_point`

For degenerate edges (zero length), returns a zero normal.
"""
function edge_outward_normal(g::AbstractManifoldMesh, edge_id::Int, cell_id::Int)
    error("$(typeof(g)) must implement `edge_outward_normal`")
end

# -- Boundary Markers --

"""
    boundary_nodes(g, marker) -> Vector{Int}

Node IDs on the boundary with the given `marker`. Returns empty for closed manifolds (e.g., full sphere).
"""
function boundary_nodes(g::AbstractManifoldMesh, marker)
    error("$(typeof(g)) must implement `boundary_nodes`")
end

"""
    boundary_edges(g, marker) -> Vector{Int}

Edge IDs on the boundary with the given `marker`. Returns empty for closed manifolds (e.g., full sphere).
"""
function boundary_edges(g::AbstractManifoldMesh, marker)
    error("$(typeof(g)) must implement `boundary_edges`")
end

# -- Trait Functions --

"""
    CellTypeStyle(g) -> CellTypeStyle
    CellTypeStyle(::Type{MyMesh}) -> CellTypeStyle

Returns the cell type trait for mesh `g` (or mesh type):
- `IsUniform{K}`: all cells have exactly K nodes
- `IsMixed{MAX_K}`: cells have varying node counts, capped at MAX_K
"""
CellTypeStyle(g::AbstractManifoldMesh) = CellTypeStyle(typeof(g))

"""
    PatchStyle(g) -> PatchStyle
    PatchStyle(::Type{MyMesh}) -> PatchStyle

Returns the patch structure trait for mesh `g` (or mesh type):
- `NoPatch`: single contiguous domain
- `MultiPatch{N}`: N independent structured patches (e.g., cubed-sphere faces)
"""
PatchStyle(g::AbstractManifoldMesh) = PatchStyle(typeof(g))

# -- Dual Mesh Functions --

"""
    dual(g) -> AbstractManifoldMesh

Returns the dual mesh of `g`, computing it lazily on first call.
Throws an error if the mesh type does not support dual construction.
"""
function dual(g::AbstractManifoldMesh)
    error("$(typeof(g)) must implement `dual` or declare no dual support")
end

# -- Patch Query Functions (for MultiPatch meshes) --

"""
    cell_face(g, cell_id) -> Int

Returns the face/patch index (1-based) that contains `cell_id`.
Only defined for meshes with `PatchStyle(g) == MultiPatch{N}`.
"""
function cell_face(g::AbstractManifoldMesh, cell_id::Int)
    error("$(typeof(g)) does not support patch queries")
end

"""
    cell_local_2d(g, cell_id) -> Tuple{Int, Int}

Returns the local 2D indices `(i, j)` of `cell_id` within its face/patch.
Only defined for meshes with `PatchStyle(g) == MultiPatch{N}`.
"""
function cell_local_2d(g::AbstractManifoldMesh, cell_id::Int)
    error("$(typeof(g)) does not support patch queries")
end
