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
using Manifolds   # Sphere

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
    _twin_nodes::Dict{Int, Vector{Int}}   # node id -> ids at the same position
    _twin_edges::Dict{Int, Vector{Int}}   # edge id -> edges with the same endpoints
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

    # Geometric twins: periodic/foreign tables may represent one geometric
    # vertex or edge by several node ids (poles, lon=0/360 seams). Location
    # tie-breaking must union incident cells across twins. Positions are
    # matched on a 1e-12 grid of the unit sphere — duplicated features in
    # real files are bit-identical points.
    _twin_nodes = Dict{Int, Vector{Int}}()
    pos_groups = Dict{Tuple{Int, Int, Int}, Vector{Int}}()
    for n in 1:n_nodes
        u = normalize(points[n])
        key = (round(Int, u[1] * 1e12), round(Int, u[2] * 1e12),
            round(Int, u[3] * 1e12))
        push!(get!(pos_groups, key, Int[]), n)
    end
    for group in values(pos_groups)
        length(group) > 1 && for n in group
            _twin_nodes[n] = group
        end
    end
    _twin_edges = Dict{Int, Vector{Int}}()
    edge_groups = Dict{Tuple{Tuple{Int, Int, Int}, Tuple{Int, Int, Int}}, Vector{Int}}()
    for e in 1:(topo.n_edges)
        i, j = getindex_fixed(topo.edge_nodes, e, Val(2))
        ui, uj = normalize(points[i]), normalize(points[j])
        ki = (round(Int, ui[1] * 1e12), round(Int, ui[2] * 1e12),
            round(Int, ui[3] * 1e12))
        kj = (round(Int, uj[1] * 1e12), round(Int, uj[2] * 1e12),
            round(Int, uj[3] * 1e12))
        key = minmax(ki, kj)
        push!(get!(edge_groups, key, Int[]), e)
    end
    for group in values(edge_groups)
        length(group) > 1 && for e in group
            _twin_edges[e] = group
        end
    end

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
        unit_corners = [normalize(p) for p in v]
        cbar = Manifolds.mean(unit, unit_corners)
        cunit = normalize(SVector{3, Float64}(cbar))
        # Hemisphere-spanning cells make the Riemannian mean collapse to the
        # zero vector (normalize then yields NaN) and l'Huilier produce
        # garbage areas; without a guard the poisoned caches only surface
        # later, far from the cause, in the locate k-d tree. The symmetric
        # corner sum is ~0 exactly when no open hemisphere contains all
        # corners (antipodal corners); coincident corners are fine (sum = K).
        (isfinite(area) && all(isfinite, cunit) &&
         norm(sum(unit_corners)) > 1e-8) || throw(ArgumentError(
            "cell $c has degenerate geometry (non-finite area or centroid, or corners spanning a hemisphere); cells must be convex spherical polygons within a hemisphere"))
        _cell_centroids[c] = R * P(cunit)
    end

    return UnstructuredMesh{typeof(M), P, MAX_K}(M, R, points, _cell_nodes,
        _cell_volumes, _cell_centroids, topo.edge_nodes, topo.edge_cells,
        topo.node_cells, topo.cell_cells, topo.node_edges, topo.cell_edges,
        _twin_nodes, _twin_edges,
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

# -- Point location --
#
# Two-phase: (1) k-d tree over cell centroids (3D Euclidean — same coarse
# semantics as UXarray's ball tree), (2) exact spherical
# point-in-convex-polygon. The k-d tree lives on representation coordinates,
# which is valid for sphere-like embeddings; periodic or curved
# representations would need a per-manifold index (override door left open).

function _locate_kdtree(g::UnstructuredMesh)
    if g._locate_index[] === nothing
        cents = reduce(hcat, g._cell_centroids)   # 3 x num_cells
        g._locate_index[] = NearestNeighbors.KDTree(cents)
    end
    return g._locate_index[]
end

"""Sign of the side of `p` w.r.t. the great circle through unit vectors `n1, n2`.
Values within 1e-12 of the plane count as on it (0.0) — FP robustness for
boundary points reconstructed through lat/lon round-trips."""
@inline function _gc_side(n1::SVector{3, Float64}, n2::SVector{3, Float64},
        p::SVector{3, Float64})
    d = dot(cross(n1, n2), p)
    return abs(d) < 1e-12 ? 0.0 : sign(d)
end

function _cell_unit_corners(g::UnstructuredMesh, cell_id::Int)
    ns = cell_nodes(g, cell_id)
    K = length(ns)
    return [normalize(SVector{3, Float64}(g.nodes[ns[k]])) for k in 1:K]
end

"""
    _point_in_convex_poly(g, cell_id, q) -> Bool

Exact spherical point-in-convex-polygon test for a cell of any arity. The
orientation reference is the cell's own centroid side: for a convex cyclic
polygon, `dot(cross(v_k, v_{k+1}), centroid)` has one sign for every edge; `q`
is inside when each edge side matches it (0 counts as on-edge/inside).
Degenerate edges (coincident corners, as in polar caps) impose no constraint
and never fix the orientation.
"""
function _point_in_convex_poly(g::UnstructuredMesh, cell_id::Int, q::SVector{3, Float64})
    v = _cell_unit_corners(g, cell_id)
    K = length(v)
    c = normalize(SVector{3, Float64}(g._cell_centroids[cell_id]))
    orient = 0.0
    for k in 1:K
        n = cross(v[k], v[mod1(k + 1, K)])
        norm(n) < 1e-12 && continue     # degenerate edge (coincident corners,
        # e.g. polar caps): no half-space constraint
        if orient == 0.0
            orient = sign(dot(n, c))    # first non-degenerate edge fixes orientation
        end
        s = _gc_side(v[k], v[mod1(k + 1, K)], q)
        (s == orient || s == 0) || return false
    end
    return true
end

"""
    _incident_containers(g, cell_id, q) -> Int

Smallest cell id among the cells containing `q`, given `q` lies on the
boundary of `cell_id`: unions the cells across every edge and vertex of
`cell_id` that `q` touches — including geometric twins (seam/pole duplicates
of the same feature carry distinct node/edge ids) — then takes the minimum
container. Exact under exact arithmetic; vertices are matched with a
1e-12 dot-product tolerance.
"""
function _incident_containers(g::UnstructuredMesh, cell_id::Int, q::SVector{3, Float64})
    v = _cell_unit_corners(g, cell_id)
    K = length(v)
    ns = collect(cell_nodes(g, cell_id))
    cands = Set{Int}([cell_id])
    for k in 1:K
        n1, n2 = v[k], v[mod1(k + 1, K)]
        if _gc_side(n1, n2, q) == 0               # q on this edge's great circle
            e = g._cell_edges.values[g._cell_edges.offsets[cell_id] - 1 + k]
            for e2 in get(g._twin_edges, e, (e,))
                for c2 in g._edge_cells[e2]
                    push!(cands, c2)
                end
            end
        end
        if dot(n1, q) > 1 - 1e-12                 # q coincides with vertex k
            for n2 in get(g._twin_nodes, ns[k], (ns[k],))
                for c2 in g._node_cells[n2]
                    push!(cands, c2)
                end
            end
        end
    end
    return minimum(c for c in cands if _point_in_convex_poly(g, c, q))
end

@inline function _locate_cell(g::UnstructuredMesh, lat::Real, lon::Real)
    -90 <= lat <= 90 || throw(ArgumentError("lat $lat out of [-90, 90]"))
    q = _latlon_to_cartesian(lat, Float64(lon), 1.0)
    tree = _locate_kdtree(g)
    n = num_cells(g)
    k = 1
    while k <= n
        idxs, _ = NearestNeighbors.knn(tree, q, k, true)
        for cid in idxs
            if _point_in_convex_poly(g, cid, q)
                # strictly interior? then cid is the answer; otherwise resolve
                # shared-boundary ties to the smallest incident container
                v = _cell_unit_corners(g, cid)
                K = length(v)
                onboundary = any(k2 -> _gc_side(v[k2], v[mod1(k2 + 1, K)], q) == 0, 1:K) ||
                             any(k2 -> dot(v[k2], q) > 1 - 1e-12, 1:K)
                return onboundary ? _incident_containers(g, cid, q) : cid
            end
        end
        k = k == n ? n + 1 : min(2 * k, n)
    end
    throw(ArgumentError("point (lat=$lat, lon=$lon) lies outside the mesh"))
end

locate_cell(g::UnstructuredMesh, lat::Real, lon::Real) = _locate_cell(g, lat, lon)

# -- Interpolation --
#
# The shared 4-node corner solver lives in locate.jl; `_wachspress_weights`
# lives HERE because its signature mentions `UnstructuredMesh`, which does not
# exist yet when locate.jl is included (module include order).

function _cell_local_coords(g::UnstructuredMesh, cell_id::Int, lat::Real, lon::Real)
    _local_coords_via_corners(g, cell_id, lat, lon)
end

"""
    _wachspress_weights(g::UnstructuredMesh, cell_id, lat, lon) -> weights

Wachspress coordinates of the query direction within a convex K-gon cell
(K != 4 path of interpolation): 2D cross-product areas on the gnomonic
projection about the cell's symmetric corner mean (order-independent). With `A_k` the signed area of triangle
`(v_k, v_{k+1}, q)` and `C_k` the area of the vertex wedge
`(v_{k-1}, v_k, v_{k+1})`, the weights are `W_k ∝ C_k / (A_{k-1} A_k)`,
normalized to sum 1. Reduces to exact barycentric coordinates on triangles.
"""
function _wachspress_weights(g::UnstructuredMesh, cell_id::Int,
        lat::Real, lon::Real)
    # Projection center: the SYMMETRIC corner mean, NOT the cached centroid.
    # `Manifolds.mean` (GeodesicInterpolation) is order-dependent sequential
    # slerp — 5-8 deg off the symmetric center for large cells — which would
    # make weights depend on the file's corner order. Any interior center is
    # valid for the gnomonic projection; the symmetric mean is order-free.
    v = _cell_unit_corners(g, cell_id)
    K = length(v)
    c = normalize(sum(v))
    if c[1]^2 + c[2]^2 < 1e-14
        # center exactly at a pole: the standard east/north basis degenerates;
        # any orthonormal frame is exact (the constructions are frame-invariant)
        east = SVector(1.0, 0.0, 0.0)
        north = SVector(0.0, 1.0, 0.0)
    else
        north = SVector(-c[1] * c[3], -c[2] * c[3], c[1]^2 + c[2]^2)
        north = north / norm(north)
        east = SVector(-c[2], c[1], 0.0)
        east = east / norm(east)
    end
    proj(p) = (dot(p, east), dot(p, north))
    θ = π / 2 - deg2rad(Float64(lat))
    φ = deg2rad(mod(Float64(lon), 360.0))
    sθ, cθ = sincos(θ)
    sφ, cφ = sincos(φ)
    q = proj(SVector{3, Float64}(sθ * cφ, sθ * sφ, cθ))
    p = [proj(v[k]) for k in 1:K]
    tri2(a, b, o) = (b[1] - a[1]) * (o[2] - a[2]) - (b[2] - a[2]) * (o[1] - a[1])
    Cs = [tri2(p[mod1(k - 1, K)], p[k], p[mod1(k + 1, K)]) for k in 1:K]
    orient = sign(Cs[1])
    A = [max(orient * tri2(p[k], p[mod1(k + 1, K)], q), 1e-14) for k in 1:K]
    C = abs.(Cs)
    W = [C[k] / (A[mod1(k - 1, K)] * A[k]) for k in 1:K]
    total = sum(W)
    return ntuple(k -> W[k] / total, K)
end

# K-aware override of the quad-bilinear default (locate.jl): quads keep the
# bilinear corner solver (consistency with the parametric grids — this is what
# the LatLonGrid oracle tests exercise); other arities use Wachspress.
function interpolation_weights(g::UnstructuredMesh, cell_id::Int,
        lat::Real, lon::Real)
    nodes = cell_nodes(g, cell_id)
    if length(nodes) == 4
        s, t = _cell_local_coords(g, cell_id, lat, lon)
        return (nodes, _bilinear_weights(s, t))
    end
    return (nodes, _wachspress_weights(g, cell_id, lat, lon))
end
