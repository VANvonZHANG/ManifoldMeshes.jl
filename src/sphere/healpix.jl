using Manifolds: Sphere
using StaticArrays: SMatrix

# -- Struct --

struct HEALPixGrid{M <: AbstractManifold} <: AbstractManifoldMesh{M}
    manifold::M
    nside::Int
    R::Float64
    ordering::Symbol
    nodes::Vector{SVector{3, Float64}}
    cell_volumes::Vector{Float64}
    cell_centroids::Vector{SVector{3, Float64}}
    _cell_nodes::Vector{NTuple{4, Int}}
    _dual::Base.RefValue{Union{Nothing, AbstractManifoldMesh{M}}}
    rotation::SMatrix{3, 3, Float64, 9}
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

# -- Constructor --

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
    _cell_nodes = Vector{NTuple{4, Int}}(undef, n_cells)

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

            _cell_nodes[cell_id] = (sw_id, se_id, ne_id, nw_id)
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

    # Scale volumes to enforce exact area conservation (compensates for
    # non-conforming gaps/overlaps in the approximate corner reconstruction).
    total_area = sum(cell_volumes)
    if total_area > 0
        scale = 4π * R^2 / total_area
        for cid in 1:n_cells
            cell_volumes[cid] *= scale
        end
    end

    return HEALPixGrid{typeof(M)}(
        M, nside, R, ordering, nodes,
        cell_volumes, cell_centroids, _cell_nodes,
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

# -- Geometry --

function node_coordinates(g::HEALPixGrid, node_id::Int)
    _check_node_id(g, node_id)
    return g.nodes[node_id]
end

function cell_volume(g::HEALPixGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g.cell_volumes[cell_id]
end

function cell_centroid(g::HEALPixGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g.cell_centroids[cell_id]
end

# -- Topology --

function cell_nodes(g::HEALPixGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g._cell_nodes[cell_id]
end

function cell_cells(g::HEALPixGrid, cell_id::Int)
    error("not yet implemented")
end

function node_cells(g::HEALPixGrid, node_id::Int)
    error("not yet implemented")
end

function cell_edges(g::HEALPixGrid, cell_id::Int)
    error("not yet implemented")
end

# -- Edge Stubs --
# TODO(phase3): implement after vertex-based reconstruction

function edge_length(g::HEALPixGrid, edge_id::Int)
    error("not yet implemented")
end

function edge_midpoint(g::HEALPixGrid, edge_id::Int)
    error("not yet implemented")
end

function edge_outward_normal(g::HEALPixGrid, edge_id::Int, cell_id::Int)
    error("not yet implemented")
end

# -- Boundary --

boundary_nodes(g::HEALPixGrid, marker) = Int[]
boundary_edges(g::HEALPixGrid, marker) = Int[]
