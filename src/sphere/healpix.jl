using Manifolds: Sphere
using StaticArrays: SMatrix

# -- Struct --

"""
    HEALPixGrid{M<:ManifoldsBase.AbstractManifold}

Hierarchical Equal Area iso-Latitude Pixelization grid.
Cells are `IsUniform{4}` quads with `IsSemiGrid` topology.
All cells have equal area; resolution is controlled by `nside` (Nside×Nside×12 cells total).
"""
struct HEALPixGrid{M <: AbstractManifold} <: AbstractManifoldMesh{M}
    manifold::M
    nside::Int
    R::Float64
    ordering::Symbol
    nodes::Vector{SVector{3, Float64}}
    cell_volumes::Vector{Float64}
    cell_centroids::Vector{SVector{3, Float64}}
    _ring_to_nested_perm::Vector{Int}   # length 12*nside²; [ipring+1] = nested_idx (0-indexed)
    _cell_nodes::CSRMapping
    _cell_edges::CSRMapping
    _cell_cells::CSRMapping
    _edge_nodes::CSRMapping
    _edge_cells::CSRMapping
    _node_edges::CSRMapping
    _node_cells::CSRMapping
    _dual::Base.RefValue{Union{Nothing, AbstractManifoldMesh{M}}}
    rotation::SMatrix{3, 3, Float64, 9}
end

# -- Internal: Ring-to-Nested Permutation --

# Standard HEALPix base pixel layout (0-indexed):
#      0  1
#   2  3  4  5
#   6  7  8  9
#     10 11

# Nested ordering: cells are ordered by base pixel, then by Morton index within each base pixel.
# Nested index = base * nside^2 + morton_index(ix, iy) (0-indexed).

"""
    _morton_decode(m::Int) -> (ix::Int, iy::Int)

Decode a Morton (Z-order) index into 2D grid coordinates.
Interleaves bits: m = ...b2a2b1a1b0a0 -> ix = ...a2a1a0, iy = ...b2b1b0.
"""
function _morton_decode(m::Int)
    ix = 0
    iy = 0
    bit = 0
    while m > 0
        ix |= (m & 1) << bit
        m >>= 1
        iy |= (m & 1) << bit
        m >>= 1
        bit += 1
    end
    return (ix, iy)
end

"""
    _morton_encode(ix::Int, iy::Int) -> m::Int

Encode 2D grid coordinates into a Morton (Z-order) index.
"""
function _morton_encode(ix::Int, iy::Int)
    m = 0
    bit = 0
    while ix > 0 || iy > 0
        m |= (ix & 1) << (2 * bit)
        m |= (iy & 1) << (2 * bit + 1)
        ix >>= 1
        iy >>= 1
        bit += 1
    end
    return m
end

"""
    _nested_cell_center(nside, base, ix, iy) -> SVector{3,Float64}

Compute the 3D Cartesian coordinates of a HEALPix cell center in nested ordering.

Parameters:
- nside: HEALPix resolution parameter
- base: base pixel index (0-11)
- ix, iy: integer coordinates within the base pixel, in [0, nside-1]

Returns the cell center as a unit vector (R=1).

This implements the standard HEALPix nested ordering geometry.
"""
function _nested_cell_center(nside::Int, base::Int, ix::Int, iy::Int)
    # Convert (base, ix, iy) to nested index then to (theta, phi)
    morton = _morton_encode(ix, iy)
    nested_idx = base * nside * nside + morton  # 0-indexed
    return _nested_to_ang(nside, nested_idx)
end

"""
    _nested_to_ang(nside, nested_idx) -> SVector{3,Float64}

Convert a 0-indexed nested cell index to 3D Cartesian coordinates on the unit sphere.
Implements the standard HEALPix nested ordering geometry.
"""
function _nested_to_ang(nside::Int, nested_idx::Int)
    ipring = _nested_to_ring(nside, nested_idx)
    return _ring_to_ang(nside, ipring)
end

"""
    _ring_to_ang(nside, ipring) -> SVector{3,Float64}

Convert a 0-indexed ring cell index to 3D Cartesian coordinates on the unit sphere.
"""
function _ring_to_ang(nside::Int, ipring::Int)
    n_cells = 12 * nside * nside
    @assert 0 <= ipring < n_cells

    # Find the ring
    # North polar cap: rings 1 to nside, ring r has 4*r cells
    # Equatorial: rings nside+1 to 3*nside, each has 4*nside cells
    # South polar cap: rings 3*nside+1 to 4*nside-1, ring r has 4*(4*nside-r) cells

    n_rings = 4 * nside - 1

    # Find ring by accumulating
    ring = 1
    remaining = ipring
    while ring <= n_rings
        if ring <= nside
            n_in_ring = 4 * ring
        elseif ring <= 3 * nside
            n_in_ring = 4 * nside
        else
            n_in_ring = 4 * (4 * nside - ring)
        end

        if remaining < n_in_ring
            break
        end
        remaining -= n_in_ring
        ring += 1
    end

    # remaining is the 0-indexed cell within the ring
    i_in_ring = remaining

    # Compute theta (colatitude)
    if ring <= nside
        # North polar cap
        cos_theta = 1.0 - ring^2 / (3.0 * nside^2)
    elseif ring <= 3 * nside
        # Equatorial belt
        cos_theta = (4.0 * nside - 2.0 * ring) / (3.0 * nside)
    else
        # South polar cap
        j = n_rings - ring + 1
        cos_theta = -(1.0 - j^2 / (3.0 * nside^2))
    end
    theta = acos(cos_theta)

    # Compute phi (longitude)
    if ring <= nside
        n_in_ring = 4 * ring
    elseif ring <= 3 * nside
        n_in_ring = 4 * nside
    else
        n_in_ring = 4 * (4 * nside - ring)
    end
    phi = 2π * (i_in_ring + 0.5) / n_in_ring

    x = sin(theta) * cos(phi)
    y = sin(theta) * sin(phi)
    z = cos(theta)

    return SVector(x, y, z)
end

# -- Internal: Bounds Checking --

@inline function _check_cell_id(g::HEALPixGrid, cell_id::Int)
    @boundscheck 1 <= cell_id <= num_cells(g) ||
                 throw(BoundsError("cell_id $cell_id out of range [1, $(num_cells(g))]"))
    nothing
end

@inline function _check_node_id(g::HEALPixGrid, node_id::Int)
    @boundscheck 1 <= node_id <= num_nodes(g) ||
                 throw(BoundsError("node_id $node_id out of range [1, $(num_nodes(g))]"))
    nothing
end

@inline function _check_edge_id(g::HEALPixGrid, edge_id::Int)
    @boundscheck 1 <= edge_id <= num_edges(g) ||
                 throw(BoundsError("edge_id $edge_id out of range [1, $(num_edges(g))]"))
    nothing
end

# -- Constructor --

"""
    HEALPixGrid(; nside::Int, ordering::Symbol = :ring, rotation = I, R::Float64 = 1.0)

Construct a HEALPix grid on the sphere.

# Arguments
- `nside`: Resolution parameter (≥ 1). Cell count = 12 × nside²
- `ordering`: `:ring` (default) or `:nested`
- `rotation`: 3×3 rotation matrix applied to all nodes
- `R`: Sphere radius (default 1.0)

HEALPix provides an equal-area hierarchical subdivision of the sphere
into 12 base pixels, each divided into nside² quadrilateral cells.
"""
function HEALPixGrid(; nside::Int, ordering::Symbol = :ring,
        rotation = SMatrix{3, 3, Float64, 9}(I), R::Float64 = 1.0)
    nside >= 1 || throw(ArgumentError("nside must be >= 1, got $nside"))
    ordering in (:ring, :nested) ||
        throw(ArgumentError("ordering must be :ring or :nested, got $ordering"))
    R > 0 || throw(ArgumentError("R must be positive, got $R"))
    rotation = convert(SMatrix{3, 3, Float64, 9}, rotation)

    M = Sphere(2)
    n_cells = 12 * nside * nside

    n_rings = 4 * nside - 1

    # --- Step 1: Generate cell centers on rings ---
    # ring_centers[ring][i] = center of cell i in ring (1-indexed)
    ring_centers = Vector{Vector{SVector{3, Float64}}}(undef, n_rings)
    ring_counts = Vector{Int}(undef, n_rings)

    for ring in 1:n_rings
        if ring <= nside
            # North polar cap
            n_in_ring = 4 * ring
            cos_theta = 1.0 - ring^2 / (3.0 * nside^2)
            theta = acos(cos_theta)
        elseif ring <= 3 * nside
            # Equatorial belt
            n_in_ring = 4 * nside
            cos_theta = (4.0 * nside - 2.0 * ring) / (3.0 * nside)
            theta = acos(cos_theta)
        else
            # South polar cap
            j = n_rings - ring + 1
            n_in_ring = 4 * j
            cos_theta = -(1.0 - j^2 / (3.0 * nside^2))
            theta = acos(cos_theta)
        end

        ring_counts[ring] = n_in_ring
        centers = SVector{3, Float64}[]

        for i in 1:n_in_ring
            phi = 2π * (i - 0.5) / n_in_ring
            x = sin(theta) * cos(phi)
            y = sin(theta) * sin(phi)
            z = cos(theta)
            p = rotation * SVector(x, y, z)
            push!(centers, SVector{3, Float64}(p))
        end
        ring_centers[ring] = centers
    end

    # --- Step 2: Compute corner vertices for each cell ---
    # For each cell, the 4 corners are at the intersections of ring boundaries
    # and sector boundaries. Each corner is the normalized average of the
    # adjacent cell centers that surround that corner.
    # For boundary corners at poles, use the pole point + adjacent centers.

    # Helper to get a ring center with periodic wrap
    function get_center(ring, i)
        n = ring_counts[ring]
        idx = mod(i - 1, n) + 1
        return ring_centers[ring][idx]
    end

    north_pole = rotation * SVector(0.0, 0.0, 1.0)
    south_pole = rotation * SVector(0.0, 0.0, -1.0)

    # Deduplicate corners using a Dict with rounded coordinates as keys
    corner_dict = Dict{Tuple{Int, Int, Int}, Int}()
    nodes = SVector{3, Float64}[]
    _cell_nodes = CSRMapping(n_cells, 4)

    function add_corner(v::SVector{3, Float64})
        # Round to 12 decimal places for deduplication (~1e-12 m precision at R=1)
        key = (round(Int, v[1] * 1e12), round(Int, v[2] * 1e12), round(Int, v[3] * 1e12))
        if haskey(corner_dict, key)
            return corner_dict[key]
        end
        idx = length(nodes) + 1
        push!(nodes, v)
        corner_dict[key] = idx
        return idx
    end

    function normalized_mean(vectors)
        s = sum(vectors)
        n = norm(s)
        if n < 1e-15
            error("normalized_mean: zero vector")
        end
        return R * (s / n)
    end

    cell_id = 1
    for ring in 1:n_rings
        n_in_ring = ring_counts[ring]
        for i in 1:n_in_ring
            c = get_center(ring, i)
            w = get_center(ring, i - 1)
            e = get_center(ring, i + 1)

            # Determine north and south neighbors
            if ring > 1
                n_ring = ring - 1
                n_count = ring_counts[n_ring]
                # Find the cell in the north ring that is closest in longitude
                # The north ring has fewer or equal cells
                # For HEALPix ring ordering, cell i in ring maps to approximately
                # the same longitudinal sector in the north ring
                ratio = n_count / n_in_ring
                i_north = clamp(round(Int, (i - 0.5) * ratio + 0.5), 1, n_count)
                n = get_center(n_ring, i_north)
                nw = get_center(n_ring, i_north - 1)
                ne = get_center(n_ring, i_north + 1)
            else
                # North pole
                n = nothing
            end

            if ring < n_rings
                s_ring = ring + 1
                s_count = ring_counts[s_ring]
                ratio = s_count / n_in_ring
                i_south = clamp(round(Int, (i - 0.5) * ratio + 0.5), 1, s_count)
                s = get_center(s_ring, i_south)
                sw_s = get_center(s_ring, i_south - 1)
                se_s = get_center(s_ring, i_south + 1)
            else
                # South pole
                s = nothing
            end

            # Compute 4 corners (SW, SE, NE, NW in local cell coordinates)
            # For ring 1 (north polar cap), the "north" side is the pole.
            # For ring n_rings (south polar cap), the "south" side is the pole.

            # SW corner: on the southern boundary, western side
            if ring == n_rings
                # South polar cap: SW corner is the south pole
                sw = south_pole
            else
                sw = normalized_mean([c, w, s, sw_s])
            end

            # SE corner: on the southern boundary, eastern side
            if ring == n_rings
                # South polar cap: SE corner is the south pole
                se = south_pole
            else
                se = normalized_mean([c, e, s, se_s])
            end

            # NW corner: on the northern boundary, western side
            if ring == 1
                # North polar cap: NW corner is the north pole
                nw_corner = north_pole
            else
                nw_corner = normalized_mean([c, w, n, nw])
            end

            # NE corner: on the northern boundary, eastern side
            if ring == 1
                # North polar cap: NE corner is the north pole
                ne_corner = north_pole
            else
                ne_corner = normalized_mean([c, e, n, ne])
            end

            sw_id = add_corner(sw)
            se_id = add_corner(se)
            ne_id = add_corner(ne_corner)
            nw_id = add_corner(nw_corner)

            base = _cell_nodes.offsets[cell_id] - 1
            _cell_nodes.values[base + 1] = sw_id
            _cell_nodes.values[base + 2] = se_id
            _cell_nodes.values[base + 3] = ne_id
            _cell_nodes.values[base + 4] = nw_id
            cell_id += 1
        end
    end

    # --- Step 3: Compute cell volumes and centroids from corner vertices ---
    cell_volumes = Vector{Float64}(undef, n_cells)
    cell_centroids = Vector{SVector{3, Float64}}(undef, n_cells)

    for cid in 1:n_cells
        cn = _cell_nodes[cid]
        A = nodes[cn[1]]   # SW
        B = nodes[cn[2]]   # SE
        C = nodes[cn[3]]   # NE
        D = nodes[cn[4]]   # NW

        area = spherical_triangle_area(R, A, B, C) +
               spherical_triangle_area(R, A, C, D)
        cell_volumes[cid] = area

        verts = [A, B, C, D]
        c = Manifolds.mean(M, verts)
        cell_centroids[cid] = SVector{3, Float64}(c)
    end

    # --- Derive edges from cell-node connectivity ---
    edge_map = Dict{Tuple{Int, Int}, Int}()
    _cell_edges = NTuple{4, Int}[]

    for cell_id in 1:n_cells
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

    n_edges = length(edge_map)
    _cell_edges_csr = CSRMapping(n_cells, 4)
    for cid in 1:n_cells
        ce = _cell_edges[cid]
        base = _cell_edges_csr.offsets[cid] - 1
        for j in 1:4
            _cell_edges_csr.values[base + j] = ce[j]
        end
    end
    _cell_edges = _cell_edges_csr

    _edge_nodes = CSRMapping(n_edges, 2)
    for ((n1, n2), edge_id) in edge_map
        base = _edge_nodes.offsets[edge_id] - 1
        _edge_nodes.values[base + 1] = n1
        _edge_nodes.values[base + 2] = n2
    end

    # --- Derive cell neighbors from edge sharing ---
    _edge_cells_tmp = [Int[] for _ in 1:n_edges]
    for cell_id in 1:n_cells
        for e in getindex_fixed(_cell_edges, cell_id, Val(4))
            push!(_edge_cells_tmp[e], cell_id)
        end
    end

    _cell_cells = CSRMapping(n_cells, 4)
    for cell_id in 1:n_cells
        ce = getindex_fixed(_cell_edges, cell_id, Val(4))
        neighbors = Int[]
        for e in ce
            adj = _edge_cells_tmp[e]
            if length(adj) == 1
                # Boundary edge or self-loop — no neighbor
                push!(neighbors, 0)
            else
                other = first(c for c in adj if c != cell_id)
                push!(neighbors, other)
            end
        end
        base = _cell_cells.offsets[cell_id] - 1
        for (j, nbr) in enumerate(neighbors)
            _cell_cells.values[base + j] = nbr
        end
    end

    # --- Derive edge → cells ---
    edge_cell_counts = fill(0, n_edges)
    for cid in 1:n_cells
        for eid in getindex_fixed(_cell_edges, cid, Val(4))
            edge_cell_counts[eid] += 1
        end
    end
    _edge_cells, ptrs = CSRMapping(n_edges, edge_cell_counts)

    for cid in 1:n_cells
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
        _node_edges.values[ptrs[n1]] = eid;
        ptrs[n1] += 1
        _node_edges.values[ptrs[n2]] = eid;
        ptrs[n2] += 1
    end

    # --- Derive node → cells ---
    node_cell_counts = fill(0, n_nodes)
    for cell_id in 1:n_cells
        for node_id in getindex_fixed(_cell_nodes, cell_id, Val(4))
            node_cell_counts[node_id] += 1
        end
    end
    _node_cells, ptrs = CSRMapping(n_nodes, node_cell_counts)

    for cell_id in 1:n_cells
        for node_id in getindex_fixed(_cell_nodes, cell_id, Val(4))
            @inbounds _node_cells.values[ptrs[node_id]] = cell_id
            @inbounds ptrs[node_id] += 1
        end
    end

    # --- Apply nested ordering permutation if requested ---
    # Build ring→nested permutation: _ring_to_nested_perm[ipring+1] = nested_idx (0-indexed)
    ring_to_nested_perm = Vector{Int}(undef, n_cells)
    for nested_idx in 0:(n_cells - 1)
        ipring = _nested_to_ring(nside, nested_idx)
        ring_to_nested_perm[ipring + 1] = nested_idx
    end

    if ordering == :nested
        # _ring_to_nested_perm[ipring+1] = nested_idx (0-indexed)
        # Build nested→ring permutation: perm[nested_id] = ring_id (1-indexed)
        perm = Vector{Int}(undef, n_cells)
        for ipring in 0:(n_cells - 1)
            nested_idx = ring_to_nested_perm[ipring + 1]
            perm[nested_idx + 1] = ipring + 1   # 1-indexed ring_id
        end

        # Build inverse permutation: inv_perm[ring_id] = nested_id (1-indexed)
        inv_perm = Vector{Int}(undef, n_cells)
        for nested_id in 1:n_cells
            inv_perm[perm[nested_id]] = nested_id
        end

        # Apply permutation to cell-related arrays
        cell_volumes = [cell_volumes[perm[i]] for i in 1:n_cells]
        cell_centroids = [cell_centroids[perm[i]] for i in 1:n_cells]

        _cell_nodes = _permute_uniform_csr(_cell_nodes, perm, Val(4))
        _cell_edges = _permute_uniform_csr(_cell_edges, perm, Val(4))

        # _cell_cells needs neighbor ID remapping via inv_perm
        nc = n_cells
        new_values = Vector{Int}(undef, nc * 4)
        for i in 1:nc
            base_old = _cell_cells.offsets[perm[i]] - 1
            base_new = (i - 1) * 4
            for j in 1:4
                old_neighbor = _cell_cells.values[base_old + j]
                new_values[base_new + j] = old_neighbor == 0 ? 0 : inv_perm[old_neighbor]
            end
        end
        _cell_cells = CSRMapping(_cell_cells.offsets[1:(nc + 1)], new_values)

        # Permute _node_cells (cell IDs change, node IDs don't)
        _node_cells_perm, _ = CSRMapping(n_nodes, node_cell_counts)
        for i in 1:n_nodes
            old = _node_cells[i]
            base = _node_cells_perm.offsets[i] - 1
            for j in 1:length(old)
                _node_cells_perm.values[base + j] = inv_perm[old[j]]
            end
        end
        _node_cells = _node_cells_perm
    else
        # :ring ordering — permutation not needed for locate, use empty vector
        ring_to_nested_perm = Int[]
    end

    return HEALPixGrid{typeof(M)}(
        M, nside, R, ordering, nodes,
        cell_volumes, cell_centroids, ring_to_nested_perm, _cell_nodes,
        _cell_edges, _cell_cells, _edge_nodes,
        _edge_cells, _node_edges, _node_cells,
        Ref{Union{Nothing, AbstractManifoldMesh{typeof(M)}}}(nothing),
        rotation)
end

# -- Trait Implementations --

TopologyStyle(::Type{<:HEALPixGrid}) = IsSemiGrid()
CellTypeStyle(::Type{<:HEALPixGrid}) = IsUniform{4}()
PatchStyle(::Type{<:HEALPixGrid}) = NoPatch()

has_dual(g::HEALPixGrid) = g._dual[] !== nothing

# -- Global Information --

manifold(g::HEALPixGrid) = g.manifold
num_cells(g::HEALPixGrid) = 12 * g.nside * g.nside
num_nodes(g::HEALPixGrid) = length(g.nodes)
num_edges(g::HEALPixGrid) = length(g._edge_nodes)

# -- Geometry --

function node_coordinates(g::HEALPixGrid, node_id::Int)
    _check_node_id(g, node_id)
    return g.nodes[node_id]
end

function cell_volume(g::HEALPixGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g.cell_volumes[cell_id]
end

all_cell_volumes(g::HEALPixGrid) = g.cell_volumes
all_node_coordinates(g::HEALPixGrid) = g.nodes

function cell_centroid(g::HEALPixGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g.cell_centroids[cell_id]
end

# -- Topology --

function cell_nodes(g::HEALPixGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return getindex_fixed(g._cell_nodes, cell_id, Val(4))
end

function cell_cells(g::HEALPixGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return getindex_fixed(g._cell_cells, cell_id, Val(4))
end

function node_cells(g::HEALPixGrid, node_id::Int)
    _check_node_id(g, node_id)
    return g._node_cells[node_id]
end

function cell_edges(g::HEALPixGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return getindex_fixed(g._cell_edges, cell_id, Val(4))
end

function edge_nodes(g::HEALPixGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    return getindex_fixed(g._edge_nodes, edge_id, Val(2))
end

function edge_cells(g::HEALPixGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    return g._edge_cells[edge_id]
end

function node_edges(g::HEALPixGrid, node_id::Int)
    _check_node_id(g, node_id)
    return g._node_edges[node_id]
end

# -- Edge Geometry --

function _edge_endpoints(g::HEALPixGrid, edge_id::Int)
    n1, n2 = getindex_fixed(g._edge_nodes, edge_id, Val(2))
    return (node_coordinates(g, n1), node_coordinates(g, n2))
end

function edge_length(g::HEALPixGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    n1, n2 = _edge_endpoints(g, edge_id)
    return Manifolds.distance(g.manifold, n1, n2)
end

function edge_midpoint(g::HEALPixGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    n1, n2 = _edge_endpoints(g, edge_id)
    if Manifolds.distance(g.manifold, n1, n2) < 1e-14
        return SVector{3, Float64}(n1)
    end
    return SVector{3, Float64}(Manifolds.mid_point(g.manifold, n1, n2))
end

function edge_outward_normal(g::HEALPixGrid, edge_id::Int, cell_id::Int)
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

boundary_nodes(g::HEALPixGrid, marker) = Int[]
boundary_edges(g::HEALPixGrid, marker) = Int[]

# -- Point location --

"""
    _nested_to_ring(nside, nested_idx) -> Int

Convert a 0-indexed nested cell index to a 0-indexed ring cell index.
Extracted from the conversion logic in `_nested_to_ang`.
"""
function _nested_to_ring(nside::Int, nested_idx::Int)
    n_cells = 12 * nside * nside
    @assert 0 <= nested_idx < n_cells

    npface = nside * nside
    face = div(nested_idx, npface)  # base pixel (0-11)
    ipf = mod(nested_idx, npface)   # index within face

    ix, iy = _morton_decode(ipf)

    jrll_arr = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]
    jpll_arr = [1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7]

    jr = jrll_arr[face + 1] * nside - ix - iy - 1
    nl4 = 4 * nside
    if jr < nside
        nr = jr
        n_before = 2 * nr * (nr - 1)
        kshift = 0
    elseif jr > 3 * nside
        nr = nl4 - jr
        n_before = n_cells - 2 * (nr + 1) * nr
        kshift = 0
    else
        nr = nside
        n_before = 2 * nside * (nside - 1) + (jr - nside) * nl4
        kshift = mod(jr - nside, 2)
    end

    jp = div(jpll_arr[face + 1] * nr + ix - iy + 1 + kshift, 2)
    if jp > nl4
        jp -= nl4
    end
    if jp < 1
        jp += nl4
    end

    return n_before + jp - 1  # 0-indexed ring cell ID
end

"""
    _ang2pix_ring(nside, theta, phi) -> Int

0-indexed ring pixel index containing direction (theta, phi).
Inverse of `_ring_to_ang`, matching this library's forward map convention
(Górski et al. 2005, adapted for uniform `phi = 2π(i+0.5)/n` offset).
`theta` is colatitude, `phi` longitude, both radians.

Key difference from the canonical C `ang2pix_ring`: our `_ring_to_ang` uses
`phi = 2π(i+0.5)/n_in_ring` for ALL rings (no equatorial ring-parity shift),
and classifies ring `nside` as north polar (not equatorial). The zone
threshold is therefore `z_nbound = (2*nside-1)/(3*nside)`, not `2/3`.
"""
function _ang2pix_ring(nside::Int, theta::Float64, phi::Float64)
    z = cos(theta)
    phi = mod(phi, 2π)
    z_nbound = (2 * nside - 1) / (3 * nside)
    npix = 12 * nside * nside

    if z > z_nbound
        # North polar cap: rings 1..nside (1-indexed)
        # z_r = 1 - r²/(3nside²)  →  r = nside*sqrt(3*(1-z))
        r = floor(Int, nside * sqrt(3.0 * (1.0 - z)) + 0.5)
        r = clamp(r, 1, nside)
        n_in_ring = 4 * r
        i = mod(floor(Int, phi / (2π) * n_in_ring), n_in_ring)
        n_before = 2 * r * (r - 1)
        return n_before + i
    elseif z < -z_nbound
        # South polar cap: j = 1..nside (1-indexed from south pole)
        j = floor(Int, nside * sqrt(3.0 * (1.0 + z)) + 0.5)
        j = clamp(j, 1, nside)
        n_in_ring = 4 * j
        i = mod(floor(Int, phi / (2π) * n_in_ring), n_in_ring)
        n_before = npix - 2 * j * (j + 1)
        return n_before + i
    else
        # Equatorial belt: rings nside+1..3*nside (1-indexed absolute)
        # z_r = (4nside - 2r)/(3nside)  →  r = nside*(4-3z)/2
        r = floor(Int, nside * (4.0 - 3.0 * z) / 2.0 + 0.5)
        r = clamp(r, nside + 1, 3 * nside)
        n_in_ring = 4 * nside
        i = mod(floor(Int, phi / (2π) * n_in_ring), n_in_ring)
        ncap = 2 * nside * (nside + 1)   # north cap has nside rings
        n_before = ncap + (r - nside - 1) * n_in_ring
        return n_before + i
    end
end

@inline function _locate_cell(g::HEALPixGrid, lat::Real, lon::Real)
    -90 <= lat <= 90 || throw(ArgumentError("lat $lat out of [-90, 90]"))
    theta = π / 2 - deg2rad(Float64(lat))     # colatitude
    phi = deg2rad(mod(Float64(lon), 360.0))
    ipring = _ang2pix_ring(g.nside, theta, phi)   # 0-indexed
    if g.ordering === :ring
        return ipring + 1
    else  # :nested — O(1) lookup via cached permutation
        return g._ring_to_nested_perm[ipring + 1] + 1
    end
end

"""
    _cell_local_coords(g::HEALPixGrid, cell_id, lat, lon) -> (s, t)

Local coordinates (s, t) in [0,1]² of (lat, lon) within `cell_id`.
HEALPix cells are not axis-aligned in (theta, phi), so a planar bilinear
solve on the 4 corner nodes is used (gnomonic projection + 2D Newton).
"""
function _cell_local_coords(g::HEALPixGrid, cell_id::Int, lat::Real, lon::Real)
    return _local_coords_via_corners(g, cell_id, lat, lon)
end

locate_cell(g::HEALPixGrid, lat::Real, lon::Real) = _locate_cell(g, lat, lon)
