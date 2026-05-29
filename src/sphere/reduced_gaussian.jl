using Manifolds: Sphere

# -- Internal Helpers --

"""
    _gaussian_latitudes(n::Int) -> Vector{Float64}

Compute interior node latitudes for `n` cell bands.
The `n - 1` interior latitude circles are at equal-area boundaries:
`sin(θⱼ) = -1 + 2j/n` for j = 1, ..., n-1.
Poles (±π/2) are handled separately in the constructor.
Returns latitudes from south to north.
"""
function _gaussian_latitudes(n::Int)
    lat_pts = Float64[]
    for j in 1:(n - 1)
        sin_lat = -1.0 + 2.0 * j / n
        push!(lat_pts, asin(sin_lat))
    end
    return lat_pts
end

"""
    _octahedral_lon_counts(nlat) -> Vector{Int}

Compute octahedral longitude counts for `nlat` cell bands.
Returns `nlat + 1` counts, one per node latitude circle (including poles),
ordered from south pole to north pole.
For interior circle j (1-indexed from south pole, excluding poles):
- j ≤ (nlat-1)÷2: count = max(4, 4j)
- j > (nlat-1)÷2: count = max(4, 4*(nlat - 1 - j + 1))
Poles always have 1 node each.
"""
function _octahedral_lon_counts(nlat::Int)
    counts = Int[]
    # South pole: 1 node
    push!(counts, 1)
    # Interior circles: octahedral rule (increasing from pole to equator, then decreasing)
    for j in 1:(nlat - 1)
        if j <= (nlat - 1) ÷ 2
            push!(counts, max(4, 4j))
        else
            push!(counts, max(4, 4 * (nlat - 1 - j + 1)))
        end
    end
    # North pole: 1 node
    push!(counts, 1)
    return counts
end

# -- Struct --

"""
    ReducedGaussianGrid{M<:ManifoldsBase.AbstractManifold}

Reduced Gaussian grid with Gaussian latitude bands and reduced longitude counts.
Cells are `IsUniform{4}` quads with `IsSemiGrid` topology.
Longitude count per band decreases toward the poles for near-equal cell areas.
"""
struct ReducedGaussianGrid{M <: AbstractManifold} <: AbstractManifoldMesh{M}
    manifold::M
    nlat::Int
    R::Float64
    lat_points::Vector{Float64}     # Gaussian latitudes in radians
    lon_counts::Vector{Int}         # Number of longitude cells per latitude band
    nodes::Vector{SVector{3, Float64}}
    cell_volumes::Vector{Float64}
    cell_centroids::Vector{SVector{3, Float64}}
    _cell_nodes::CSRMapping
    _cell_edges::CSRMapping
    _edge_nodes::CSRMapping
    _cell_cells::CSRMapping
    _edge_cells::CSRMapping
    _node_edges::CSRMapping
    _dual::Base.RefValue{Union{Nothing, AbstractManifoldMesh{M}}}
end

# -- Internal: Bounds Checking --

@inline function _check_cell_id(g::ReducedGaussianGrid, cell_id::Int)
    @boundscheck 1 <= cell_id <= num_cells(g) ||
                 throw(BoundsError("cell_id $cell_id out of range [1, $(num_cells(g))]"))
    nothing
end

@inline function _check_node_id(g::ReducedGaussianGrid, node_id::Int)
    @boundscheck 1 <= node_id <= num_nodes(g) ||
                 throw(BoundsError("node_id $node_id out of range [1, $(num_nodes(g))]"))
    nothing
end

@inline function _check_edge_id(g::ReducedGaussianGrid, edge_id::Int)
    @boundscheck 1 <= edge_id <= num_edges(g) ||
                 throw(BoundsError("edge_id $edge_id out of range [1, $(num_edges(g))]"))
    nothing
end

# -- Constructor --

"""
    ReducedGaussianGrid(; nlat::Int, R::Float64 = 1.0)

Construct a reduced Gaussian grid on the sphere.

# Arguments
- `nlat`: Number of latitude cell bands (≥ 2)
- `R`: Sphere radius (default 1.0)

The grid uses equal-area latitude bands with octahedral longitude reduction.
Poles are included. Cell count varies per band.
"""
function ReducedGaussianGrid(; nlat::Int, R::Float64 = 1.0)
    nlat >= 2 || throw(ArgumentError("nlat must be >= 2, got $nlat"))
    R > 0 || throw(ArgumentError("R must be positive, got $R"))

    M = Sphere(2)

    # Gaussian latitudes: nlat-1 interior boundaries + 2 poles
    interior_lats = _gaussian_latitudes(nlat)
    lat_points = vcat([-π / 2], interior_lats, [π / 2])

    # Octahedral lon counts per latitude circle
    lon_counts = _octahedral_lon_counts(nlat)

    # Build nodes on each latitude circle
    nodes = SVector{3, Float64}[]
    for j in 1:(nlat + 1)
        lat = lat_points[j]
        nlon = lon_counts[j]
        θ = π / 2 - lat  # colatitude
        sinθ, cosθ = sin(θ), cos(θ)

        for i in 1:nlon
            lon = 2π * (i - 1) / nlon
            x = R * sinθ * cos(lon)
            y = R * sinθ * sin(lon)
            z = R * cosθ
            push!(nodes, SVector(x, y, z))
        end
    end

    # Build cells between adjacent latitude circles
    # Each cell is a quadrilateral (or triangle at poles): SW, SE, NE, NW
    cell_volumes = Float64[]
    cell_centroids = SVector{3, Float64}[]
    _cell_nodes = NTuple{4, Int}[]

    for j in 1:nlat
        nlon_lower = lon_counts[j]
        nlon_upper = lon_counts[j + 1]

        lower_offset = sum(lon_counts[1:(j - 1)]; init = 0)
        upper_offset = sum(lon_counts[1:j]; init = 0)

        # Number of cells in this band = max of the two node counts
        ncells = max(nlon_lower, nlon_upper)

        for k in 1:ncells
            if nlon_lower >= nlon_upper
                # Lower is finer or equal: iterate over lower cells
                i = k
                idx_A = lower_offset + i
                idx_B = lower_offset + mod(i, nlon_lower) + 1
                A = nodes[idx_A]                              # SW
                B = nodes[idx_B]                              # SE

                ratio = nlon_upper / nlon_lower
                i_upper = ceil(Int, i * ratio)
                i_upper_next = ceil(Int, (i + 1) * ratio)
                idx_D = upper_offset + i_upper
                idx_C = upper_offset + mod(i_upper_next - 1, nlon_upper) + 1
                D = nodes[idx_D]                              # NW
                C = nodes[idx_C]                              # NE
            else
                # Upper is finer: iterate over upper cells
                i_upper = k
                idx_D = upper_offset + i_upper
                idx_C = upper_offset + mod(i_upper, nlon_upper) + 1
                D = nodes[idx_D]                              # NW
                C = nodes[idx_C]                              # NE

                ratio = nlon_lower / nlon_upper
                i = ceil(Int, i_upper * ratio)
                i_next = ceil(Int, (i_upper + 1) * ratio)
                idx_A = lower_offset + i
                idx_B = lower_offset + mod(i_next - 1, nlon_lower) + 1
                A = nodes[idx_A]                              # SW
                B = nodes[idx_B]                              # SE
            end

            area = spherical_triangle_area(R, A, B, C) +
                   spherical_triangle_area(R, A, C, D)
            push!(cell_volumes, area)

            verts = [A, B, C, D]
            c = Manifolds.mean(M, verts)
            push!(cell_centroids, SVector{3, Float64}(c))
            push!(_cell_nodes, (idx_A, idx_B, idx_C, idx_D))
        end
    end

    # --- Convert _cell_nodes to CSR ---
    num_cells = length(_cell_nodes)
    _cell_nodes_csr = CSRMapping(num_cells, 4)
    for cid in 1:num_cells
        cn = _cell_nodes[cid]
        base = _cell_nodes_csr.offsets[cid] - 1
        for j in 1:4
            _cell_nodes_csr.values[base + j] = cn[j]
        end
    end
    _cell_nodes = _cell_nodes_csr

    # --- Derive edges from cell-node connectivity ---
    edge_map = Dict{Tuple{Int, Int}, Int}()
    _cell_edges = NTuple{4, Int}[]

    for cell_id in 1:length(_cell_nodes)
        cn = _cell_nodes[cell_id]
        cell_edge_ids = Int[]
        for (a, b) in ((cn[1], cn[2]), (cn[2], cn[3]),
            (cn[3], cn[4]), (cn[4], cn[1]))
            key = a < b ? (a, b) : (b, a)
            edge_id = get!(edge_map, key) do
                length(edge_map) + 1
            end
            push!(cell_edge_ids, edge_id)
        end
        push!(_cell_edges, tuple(cell_edge_ids...))
    end

    # --- Convert _cell_edges to CSR ---
    _cell_edges_csr = CSRMapping(num_cells, 4)
    for cid in 1:num_cells
        ce = _cell_edges[cid]
        base = _cell_edges_csr.offsets[cid] - 1
        for j in 1:4
            _cell_edges_csr.values[base + j] = ce[j]
        end
    end
    _cell_edges = _cell_edges_csr

    n_edges = length(edge_map)
    _edge_nodes = CSRMapping(n_edges, 2)
    for ((n1, n2), edge_id) in edge_map
        base = _edge_nodes.offsets[edge_id] - 1
        _edge_nodes.values[base + 1] = n1
        _edge_nodes.values[base + 2] = n2
    end

    # --- Derive cell neighbors from edge sharing ---
    _edge_cells_tmp = [Int[] for _ in 1:n_edges]
    for cell_id in 1:num_cells
        for e in getindex_fixed(_cell_edges, cell_id, Val(4))
            push!(_edge_cells_tmp[e], cell_id)
        end
    end

    _cell_cells = CSRMapping(num_cells, 4)
    for cell_id in 1:num_cells
        ce = getindex_fixed(_cell_edges, cell_id, Val(4))
        neighbors = Int[]
        for e in ce
            adj = _edge_cells_tmp[e]
            n1, n2 = getindex_fixed(_edge_nodes, e, Val(2))
            if n1 == n2
                # Self-loop edge (zero-length, e.g. at poles): no neighbor across it
                push!(neighbors, 0)
            else
                # Non-degenerate edge: should be shared by exactly 2 cells
                others = filter(c -> c != cell_id, adj)
                # Use 0 sentinel if no other cell found (should not happen for non-self-loops)
                push!(neighbors, isempty(others) ? 0 : first(others))
            end
        end
        base = _cell_cells.offsets[cell_id] - 1
        for (j, nbr) in enumerate(neighbors)
            _cell_cells.values[base + j] = nbr
        end
    end

    # --- Derive edge → cells ---
    edge_cell_counts = fill(0, n_edges)
    for cid in 1:num_cells
        for eid in getindex_fixed(_cell_edges, cid, Val(4))
            edge_cell_counts[eid] += 1
        end
    end
    _edge_cells, ptrs = CSRMapping(n_edges, edge_cell_counts)

    for cid in 1:num_cells
        for eid in getindex_fixed(_cell_edges, cid, Val(4))
            pos = ptrs[eid]
            _edge_cells.values[pos] = cid
            ptrs[eid] += 1
        end
    end

    # --- Derive node → edges ---
    n_nodes = length(nodes)
    node_edge_counts = fill(0, n_nodes)
    for eid in 1:n_edges
        n1, n2 = getindex_fixed(_edge_nodes, eid, Val(2))
        node_edge_counts[n1] += 1
        node_edge_counts[n2] += 1
    end
    _node_edges, ptrs = CSRMapping(n_nodes, node_edge_counts)

    for eid in 1:n_edges
        n1, n2 = getindex_fixed(_edge_nodes, eid, Val(2))
        _node_edges.values[ptrs[n1]] = eid; ptrs[n1] += 1
        _node_edges.values[ptrs[n2]] = eid; ptrs[n2] += 1
    end

    return ReducedGaussianGrid{typeof(M)}(
        M, nlat, R, lat_points, lon_counts, nodes,
        cell_volumes, cell_centroids, _cell_nodes,
        _cell_edges, _edge_nodes, _cell_cells, _edge_cells, _node_edges,
        Ref{Union{Nothing, AbstractManifoldMesh{typeof(M)}}}(nothing))
end

# -- Trait Implementations --

TopologyStyle(::Type{<:ReducedGaussianGrid}) = IsSemiGrid()
CellTypeStyle(::Type{<:ReducedGaussianGrid}) = IsUniform{4}()
PatchStyle(::Type{<:ReducedGaussianGrid}) = NoPatch()

has_dual(g::ReducedGaussianGrid) = g._dual[] !== nothing

# -- Global Information --

manifold(g::ReducedGaussianGrid) = g.manifold
num_cells(g::ReducedGaussianGrid) = length(g.cell_volumes)
num_nodes(g::ReducedGaussianGrid) = length(g.nodes)
num_edges(g::ReducedGaussianGrid) = length(g._edge_nodes)

# -- Geometry --

function node_coordinates(g::ReducedGaussianGrid, node_id::Int)
    _check_node_id(g, node_id)
    return g.nodes[node_id]
end

function cell_volume(g::ReducedGaussianGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g.cell_volumes[cell_id]
end

all_cell_volumes(g::ReducedGaussianGrid) = g.cell_volumes
all_node_coordinates(g::ReducedGaussianGrid) = g.nodes

function cell_centroid(g::ReducedGaussianGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g.cell_centroids[cell_id]
end

# -- Topology Stubs --

function cell_nodes(g::ReducedGaussianGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return getindex_fixed(g._cell_nodes, cell_id, Val(4))
end

function cell_cells(g::ReducedGaussianGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return getindex_fixed(g._cell_cells, cell_id, Val(4))
end

function node_cells(g::ReducedGaussianGrid, node_id::Int)
    _check_node_id(g, node_id)
    cells = Int[]
    for cell_id in 1:num_cells(g)
        if node_id in cell_nodes(g, cell_id)
            push!(cells, cell_id)
        end
    end
    return cells
end

function cell_edges(g::ReducedGaussianGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return getindex_fixed(g._cell_edges, cell_id, Val(4))
end

function edge_nodes(g::ReducedGaussianGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    return getindex_fixed(g._edge_nodes, edge_id, Val(2))
end

function edge_cells(g::ReducedGaussianGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    return g._edge_cells[edge_id]
end

function node_edges(g::ReducedGaussianGrid, node_id::Int)
    _check_node_id(g, node_id)
    return g._node_edges[node_id]
end

# -- Edge Geometry --

function _edge_endpoints(g::ReducedGaussianGrid, edge_id::Int)
    n1, n2 = getindex_fixed(g._edge_nodes, edge_id, Val(2))
    return (node_coordinates(g, n1), node_coordinates(g, n2))
end

function edge_length(g::ReducedGaussianGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    n1, n2 = _edge_endpoints(g, edge_id)
    return Manifolds.distance(g.manifold, n1, n2)
end

function edge_midpoint(g::ReducedGaussianGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    n1, n2 = _edge_endpoints(g, edge_id)
    if Manifolds.distance(g.manifold, n1, n2) < 1e-14
        return SVector{3, Float64}(n1)
    end
    return SVector{3, Float64}(Manifolds.mid_point(g.manifold, n1, n2))
end

function edge_outward_normal(g::ReducedGaussianGrid, edge_id::Int, cell_id::Int)
    _check_edge_id(g, edge_id)
    _check_cell_id(g, cell_id)
    n1, n2 = _edge_endpoints(g, edge_id)

    if Manifolds.distance(g.manifold, n1, n2) < 1e-14
        midpoint = n1
        return (base_point = SVector{3, Float64}(midpoint),
            normal = zero(SVector{3, Float64}))
    end

    midpoint = Manifolds.mid_point(g.manifold, n1, n2)
    gc_normal = cross(SVector(n1), SVector(n2))
    tangent = normalize(cross(gc_normal, SVector(midpoint)))
    cell_c = cell_centroid(g, cell_id)
    cell_side = sign(dot(gc_normal, SVector(cell_c)))
    outward = cell_side * cross(tangent, SVector(midpoint))
    outward = Manifolds.project(g.manifold, midpoint, outward)

    return (base_point = SVector{3, Float64}(midpoint),
        normal = SVector{3, Float64}(outward))
end

# -- Boundary --

boundary_nodes(g::ReducedGaussianGrid, marker) = Int[]
boundary_edges(g::ReducedGaussianGrid, marker) = Int[]
