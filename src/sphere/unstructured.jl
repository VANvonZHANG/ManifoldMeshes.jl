using Manifolds: Sphere
import NearestNeighbors

"""
    UnstructuredMesh{M, P, MAX_K} <: AbstractManifoldMesh{M}

An arbitrary spherical polygon mesh stored as *data*: node coordinates plus a
padded face-node table. This is the mesh analogue of the four parametric grid
types — they *generate* everything from constructor parameters,
`UnstructuredMesh` *derives* everything from the table.

Type parameters: `M` is the manifold (v1: `Sphere` only, enforced at
construction), `P` the node point type (`SVector{3, Float64}` on the sphere),
and `MAX_K` the padded connectivity width. Cells are `IsMixed{MAX_K}` —
arbitrary convex polygons with per-cell arity 3..MAX_K; quads (K=4) interpolate
bilinearly like the parametric grids, other arities use Wachspress
coordinates. Corners are listed in **cyclic boundary order** (starting corner
irrelevant). Point location uses a lazily built k-d tree over cell centroids
(`NearestNeighbors.jl`) plus an exact spherical point-in-polygon test.

# Constructors
```julia
UnstructuredMesh(M::AbstractManifold, points::Vector{<:SVector{3,Float64}}, face_nodes; fill_value = -1)
UnstructuredMesh(node_lon, node_lat, face_nodes; R = 1.0, start_index = 0, fill_value = -1)
```
`face_nodes` is `n_cells x MAX_K`; per-cell arity is the count of leading
entries different from `fill_value`. The second (sphere convenience)
constructor takes degrees; `start_index` shifts the connectivity table to
1-based (pass `0` for UGRID files, `1` for Julia tables).
"""
struct UnstructuredMesh{M <: AbstractManifold, P, MAX_K} <: AbstractManifoldMesh{M}
    manifold::M
    R::Float64
    nodes::Vector{P}
    _cell_nodes::CSRMapping
    _cell_volumes::Vector{Float64}
    _cell_centroids::Vector{P}
    _edge_nodes::CSRMapping
    _edge_cells::CSRMapping
    _node_cells::CSRMapping
    _cell_cells::CSRMapping
    _node_edges::CSRMapping
    _cell_edges::CSRMapping
    _locate_index::Base.RefValue{Union{Nothing, NearestNeighbors.KDTree}}
    _dual::Base.RefValue{Union{Nothing, AbstractManifoldMesh{M}}}
end

@inline function _check_cell_id(g::UnstructuredMesh, cell_id::Int)
    @boundscheck 1 <= cell_id <= num_cells(g) ||
                 throw(BoundsError("cell_id $cell_id out of range [1, $(num_cells(g))]"))
    nothing
end

@inline function _check_node_id(g::UnstructuredMesh, node_id::Int)
    @boundscheck 1 <= node_id <= num_nodes(g) ||
                 throw(BoundsError("node_id $node_id out of range [1, $(num_nodes(g))]"))
    nothing
end

@inline function _check_edge_id(g::UnstructuredMesh, edge_id::Int)
    @boundscheck 1 <= edge_id <= num_edges(g) ||
                 throw(BoundsError("edge_id $edge_id out of range [1, $(num_edges(g))]"))
    nothing
end

"""
    _prepare_connectivity(face_nodes, start_index, fill_value) -> (conn, ks)

Validate a padded connectivity table: shift active (non-fill) entries to
1-based, compute per-cell arity `ks`, and enforce that each row has at least 3
active corners and that its non-fill prefix is contiguous.
"""
function _prepare_connectivity(face_nodes::AbstractMatrix{<:Integer},
        start_index::Int, fill_value::Int)
    n_cells = size(face_nodes, 1)
    n_cells > 0 || throw(ArgumentError("face_nodes must have at least one row"))
    ncol = size(face_nodes, 2)
    ncol >= 3 || throw(ArgumentError(
        "face_nodes must have at least 3 columns, got $ncol"))
    conn = Matrix{Int}(face_nodes)
    ks = Vector{Int}(undef, n_cells)
    for c in 1:n_cells
        k = 0
        while k < ncol && conn[c, k + 1] != fill_value
            conn[c, k + 1] -= (start_index - 1)
            k += 1
        end
        k >= 3 || throw(ArgumentError(
            "face_nodes row $c has $k active corners; cells need at least 3"))
        all(==(fill_value), @view conn[c, (k + 1):end]) || throw(ArgumentError(
            "face_nodes row $c mixes fill values into the active prefix; the non-fill prefix must be contiguous"))
        ks[c] = k
    end
    return conn, ks
end

function _build_unstructured(M, R::Float64, points::Vector{P}, conn::Matrix{Int},
        ks::Vector{Int}) where {P <: SVector{3, Float64}}
    n_nodes = length(points)
    n_cells = size(conn, 1)
    MAX_K = size(conn, 2)
    n_nodes > 0 || throw(ArgumentError("mesh needs at least one node"))
    for c in 1:n_cells
        row = ntuple(k -> conn[c, k], Val(MAX_K))
        active = row[1:ks[c]]
        all(i -> 1 <= i <= n_nodes, active) || throw(ArgumentError(
            "face_nodes row $c references a node outside 1:$n_nodes"))
        length(unique(active)) == ks[c] || throw(ArgumentError(
            "face_nodes row $c repeats a node; cells must have distinct corners in cyclic boundary order"))
    end

    _cell_nodes, cptrs = CSRMapping(n_cells, ks)
    for c in 1:n_cells
        for k in 1:ks[c]
            _cell_nodes.values[cptrs[c]] = conn[c, k]
            cptrs[c] += 1
        end
    end

    topo = _derive_mesh_topology(conn, ks)

    # Geometry caches: l'Huilier areas by fan triangulation from corner 1
    # (for quads this is exactly the existing A-C diagonal split, with the
    # same degenerate B-D fallback as LatLonGrid); centroids as the
    # Riemannian mean of the unit-normalized corners, rescaled by R.
    unit = Sphere(2)
    _cell_volumes = Vector{Float64}(undef, n_cells)
    _cell_centroids = Vector{P}(undef, n_cells)
    for c in 1:n_cells
        K = ks[c]
        v = [points[conn[c, k]] for k in 1:K]
        area = 0.0
        for k in 2:(K - 1)
            area += spherical_triangle_area(R, v[1], v[k], v[k + 1])
        end
        if area == 0.0 && K == 4
            area = spherical_triangle_area(R, v[2], v[3], v[4]) +
                   spherical_triangle_area(R, v[1], v[2], v[4])
        end
        _cell_volumes[c] = area
        cbar = Manifolds.mean(unit, [normalize(p) for p in v])
        _cell_centroids[c] = R * P(normalize(SVector{3, Float64}(cbar)))
    end

    return UnstructuredMesh{typeof(M), P, MAX_K}(M, R, points, _cell_nodes,
        _cell_volumes, _cell_centroids, topo.edge_nodes, topo.edge_cells,
        topo.node_cells, topo.cell_cells, topo.node_edges, topo.cell_edges,
        Ref{Union{Nothing, NearestNeighbors.KDTree}}(nothing),
        Ref{Union{Nothing, AbstractManifoldMesh{typeof(M)}}}(nothing))
end

"""
    UnstructuredMesh(M, points, face_nodes; fill_value = -1)

General entry point (v1: `M` must be a `Sphere`, checked at construction).
`R` is inferred from the mean node norm.
"""
function UnstructuredMesh(M::AbstractManifold, points::Vector{<:SVector{3, Float64}},
        face_nodes::AbstractMatrix{<:Integer}; fill_value::Integer = -1)
    M isa Sphere || throw(ArgumentError(
        "UnstructuredMesh v1 supports Sphere manifolds only, got $(typeof(M))"))
    R = sum(norm, points) / length(points)
    conn, ks = _prepare_connectivity(face_nodes, 1, Int(fill_value))
    return _build_unstructured(M, Float64(R), points, conn, ks)
end

"""
    UnstructuredMesh(node_lon, node_lat, face_nodes; R = 1.0, start_index = 0, fill_value = -1)

Sphere convenience constructor from degrees. `start_index` is the connectivity
base (0 for UGRID files, 1 for Julia tables); `fill_value` pads inactive
connectivity slots (leading non-fill prefix is the active cell).
"""
function UnstructuredMesh(node_lon::AbstractVector{<:Real},
        node_lat::AbstractVector{<:Real}, face_nodes::AbstractMatrix{<:Integer};
        R::Real = 1.0, start_index::Integer = 0, fill_value::Integer = -1)
    length(node_lon) == length(node_lat) || throw(ArgumentError(
        "node_lon has $(length(node_lon)) entries but node_lat has $(length(node_lat))"))
    points = [_latlon_to_cartesian(node_lat[i], node_lon[i], Float64(R))
              for i in eachindex(node_lon)]
    conn, ks = _prepare_connectivity(face_nodes, Int(start_index), Int(fill_value))
    return _build_unstructured(Sphere(2), Float64(R), points, conn, ks)
end

# -- Traits --

TopologyStyle(::Type{<:UnstructuredMesh}) = IsMesh()
function CellTypeStyle(::Type{<:UnstructuredMesh{M, P, MAX_K}}) where {M, P, MAX_K}
    IsMixed{MAX_K}()
end
PatchStyle(::Type{<:UnstructuredMesh}) = NoPatch()

# -- Global information and topology reads --

manifold(g::UnstructuredMesh) = g.manifold
num_cells(g::UnstructuredMesh) = length(g._cell_volumes)
num_nodes(g::UnstructuredMesh) = length(g.nodes)
num_edges(g::UnstructuredMesh) = length(g._edge_nodes)

function node_coordinates(g::UnstructuredMesh, node_id::Int)
    _check_node_id(g, node_id)
    return g.nodes[node_id]
end

function cell_nodes(g::UnstructuredMesh{M, P, MAX_K},
        cell_id::Int) where {M, P, MAX_K}
    _check_cell_id(g, cell_id)
    row = g._cell_nodes[cell_id]
    K = length(row)
    return MixedCellTopology{MAX_K}(ntuple(k -> k <= K ? row[k] : 0, Val(MAX_K)), K)
end

function cell_cells(g::UnstructuredMesh, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g._cell_cells[cell_id]
end

function node_cells(g::UnstructuredMesh, node_id::Int)
    _check_node_id(g, node_id)
    return g._node_cells[node_id]
end

function cell_edges(g::UnstructuredMesh, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g._cell_edges[cell_id]
end

function edge_nodes(g::UnstructuredMesh, edge_id::Int)
    _check_edge_id(g, edge_id)
    return getindex_fixed(g._edge_nodes, edge_id, Val(2))
end

function edge_cells(g::UnstructuredMesh, edge_id::Int)
    _check_edge_id(g, edge_id)
    return g._edge_cells[edge_id]
end

function node_edges(g::UnstructuredMesh, node_id::Int)
    _check_node_id(g, node_id)
    return g._node_edges[node_id]
end

# -- Boundary markers (closed sphere) --

boundary_nodes(g::UnstructuredMesh, marker) = Int[]
boundary_edges(g::UnstructuredMesh, marker) = Int[]

# -- Geometry --

function cell_volume(g::UnstructuredMesh, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g._cell_volumes[cell_id]
end

all_cell_volumes(g::UnstructuredMesh) = g._cell_volumes
all_cell_centroids(g::UnstructuredMesh) = g._cell_centroids
all_node_coordinates(g::UnstructuredMesh) = g.nodes

function cell_centroid(g::UnstructuredMesh, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g._cell_centroids[cell_id]
end

function _edge_endpoints(g::UnstructuredMesh, edge_id::Int)
    n1, n2 = getindex_fixed(g._edge_nodes, edge_id, Val(2))
    return (g.nodes[n1], g.nodes[n2])
end

function edge_length(g::UnstructuredMesh, edge_id::Int)
    _check_edge_id(g, edge_id)
    n1, n2 = _edge_endpoints(g, edge_id)
    return Manifolds.distance(Sphere(2), normalize(n1), normalize(n2)) * g.R
end

function edge_midpoint(g::UnstructuredMesh, edge_id::Int)
    _check_edge_id(g, edge_id)
    n1, n2 = _edge_endpoints(g, edge_id)
    u1, u2 = normalize(n1), normalize(n2)
    if Manifolds.distance(Sphere(2), u1, u2) < 1e-14
        return SVector{3, Float64}(n1)
    end
    return g.R * SVector{3, Float64}(Manifolds.mid_point(Sphere(2), u1, u2))
end

function edge_outward_normal(g::UnstructuredMesh, edge_id::Int, cell_id::Int)
    _check_edge_id(g, edge_id)
    _check_cell_id(g, cell_id)
    n1, n2 = _edge_endpoints(g, edge_id)
    u1 = normalize(SVector{3, Float64}(n1))
    u2 = normalize(SVector{3, Float64}(n2))

    # Guard: degenerate edge with coincident endpoints
    if Manifolds.distance(Sphere(2), u1, u2) < 1e-14
        return (base_point = SVector{3, Float64}(g.R * u1),
            normal = zero(SVector{3, Float64}))
    end

    m = Manifolds.mid_point(Sphere(2), u1, u2)     # unit midpoint
    gc_normal = cross(u1, u2)                      # great circle plane normal
    tangent = normalize(cross(gc_normal, m))
    c = normalize(SVector{3, Float64}(g._cell_centroids[cell_id]))
    cell_side = sign(dot(gc_normal, c))

    outward = cell_side * cross(tangent, m)
    outward = Manifolds.project(Sphere(2), m, outward)

    return (base_point = g.R * SVector{3, Float64}(m),
        normal = SVector{3, Float64}(outward))
end
