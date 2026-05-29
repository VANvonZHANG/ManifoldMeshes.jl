using Manifolds: Sphere

"""
    LatLonGrid{M<:ManifoldsBase.AbstractManifold}

Structured latitude-longitude grid on the sphere.
Cells are uniform quads with `IsGrid` topology.
Polar cells may be degenerate (poles are nodes, not cells).
"""
struct LatLonGrid{M <: AbstractManifold} <: AbstractManifoldMesh{M}
    manifold::M
    lat_edges::Vector{Float64}
    lon_edges::Vector{Float64}
    R::Float64
    nlat::Int
    nlon::Int
    nodes::Matrix{SVector{3, Float64}}
    cell_volumes::Matrix{Float64}
    cell_centroids::Matrix{SVector{3, Float64}}
    _cell_nodes::CSRMapping
    _cell_edges::CSRMapping
    _cell_cells::CSRMapping
    _edge_nodes::CSRMapping
    _edge_cells::CSRMapping
    _node_edges::CSRMapping
    _dual::Base.RefValue{Union{Nothing, AbstractManifoldMesh{M}}}
end

# -- Internal: Bounds Checking --

@inline function _check_cell_id(g::LatLonGrid, cell_id::Int)
    @boundscheck 1 <= cell_id <= num_cells(g) ||
                 throw(BoundsError("cell_id $cell_id out of range [1, $(num_cells(g))]"))
    nothing
end

@inline function _check_node_id(g::LatLonGrid, node_id::Int)
    @boundscheck 1 <= node_id <= num_nodes(g) ||
                 throw(BoundsError("node_id $node_id out of range [1, $(num_nodes(g))]"))
    nothing
end

@inline function _check_edge_id(g::LatLonGrid, edge_id::Int)
    @boundscheck 1 <= edge_id <= num_edges(g) ||
                 throw(BoundsError("edge_id $edge_id out of range [1, $(num_edges(g))]"))
    nothing
end

function LatLonGrid(; lat_edges::Vector{Float64}, lon_edges::Vector{Float64}, R::Float64 = 1.0)
    # Copy to prevent external mutation from corrupting cached state
    lat_edges = copy(lat_edges)
    lon_edges = copy(lon_edges)

    # Validate
    length(lat_edges) < 2 && throw(ArgumentError("lat_edges must have at least 2 elements"))
    length(lon_edges) < 2 && throw(ArgumentError("lon_edges must have at least 2 elements"))
    lat_edges[1] != -90.0 && throw(ArgumentError("lat_edges must start at -90"))
    lat_edges[end] != 90.0 && throw(ArgumentError("lat_edges must end at 90"))
    !issorted(lat_edges) && throw(ArgumentError("lat_edges must be ascending"))
    lon_edges[1] != 0.0 && throw(ArgumentError("lon_edges must start at 0"))
    lon_edges[end] != 360.0 && throw(ArgumentError("lon_edges must end at 360"))
    !issorted(lon_edges) && throw(ArgumentError("lon_edges must be ascending"))
    R <= 0 && throw(ArgumentError("R must be positive, got $R"))

    nlat = length(lat_edges) - 1
    nlon = length(lon_edges) - 1
    M = Sphere(2)

    # Build nodes matrix [ilat, ilon], size (nlat+1) x (nlon+1)
    nodes = Matrix{SVector{3, Float64}}(undef, nlat + 1, nlon + 1)
    for ilat in 1:(nlat + 1)
        lat = lat_edges[ilat]
        θ = deg2rad(lat)
        sinθ, cosθ = sin(θ), cos(θ)
        for ilon in 1:(nlon + 1)
            lon = lon_edges[ilon]
            φ = deg2rad(lon)
            nodes[ilat, ilon] = R * SVector(cosθ * cos(φ), cosθ * sin(φ), sinθ)
        end
    end

    # Pre-compute cell volumes via l'Huilier's formula (2 triangles per cell)
    # Split each quadrilateral into 2 triangles along the A-C diagonal.
    # If the A-C diagonal is degenerate (coincident A==C, or antipodal B==D
    # causing one triangle to span a semicircle), fall back to the B-D diagonal.
    cell_volumes = Matrix{Float64}(undef, nlat, nlon)
    for ilat in 1:nlat
        for ilon in 1:nlon
            ilon_next = ilon == nlon ? 1 : ilon + 1
            A = nodes[ilat, ilon]             # SW
            B = nodes[ilat, ilon_next]        # SE
            C = nodes[ilat + 1, ilon_next]    # NE
            D = nodes[ilat + 1, ilon]         # NW
            area_ac = spherical_triangle_area(R, A, B, C) +
                      spherical_triangle_area(R, A, C, D)
            if area_ac == 0.0
                # Degenerate A-C diagonal; use B-D diagonal instead
                area_ac = spherical_triangle_area(R, B, C, D) +
                          spherical_triangle_area(R, A, B, D)
            end
            if area_ac == 0.0
                # Both diagonals degenerate (e.g. 180° polar cells where
                # all vertices collapse to 2 antipodal points). Use lune formula.
                Δlon_rad = deg2rad(lon_edges[ilon_next] - lon_edges[ilon])
                area_ac = R^2 *
                          (sin(deg2rad(lat_edges[ilat + 1])) -
                           sin(deg2rad(lat_edges[ilat]))) * abs(Δlon_rad)
            end
            cell_volumes[ilat, ilon] = area_ac
        end
    end

    # Pre-compute cell centroids via Manifolds.mean
    cell_centroids = Matrix{SVector{3, Float64}}(undef, nlat, nlon)
    for ilat in 1:nlat
        for ilon in 1:nlon
            ilon_next = ilon == nlon ? 1 : ilon + 1
            verts = [nodes[ilat, ilon], nodes[ilat, ilon_next],
                nodes[ilat + 1, ilon_next], nodes[ilat + 1, ilon]]
            c = Manifolds.mean(M, verts)
            cell_centroids[ilat, ilon] = SVector{3, Float64}(c)
        end
    end

    n_cells = nlat * nlon
    n_nodes = (nlat + 1) * (nlon + 1)
    n_h_edges = (nlat + 1) * nlon
    n_v_edges = nlat * (nlon + 1)
    n_edges = n_h_edges + n_v_edges

    # cell → nodes (4 per cell)
    _cell_nodes = CSRMapping(n_cells, 4)
    for ilat in 1:nlat, ilon in 1:nlon
        cid = (ilat - 1) * nlon + ilon
        sw = (ilat - 1) * (nlon + 1) + ilon
        se = sw + 1
        nw = ilat * (nlon + 1) + ilon
        ne = nw + 1
        base = _cell_nodes.offsets[cid] - 1
        _cell_nodes.values[base + 1] = sw
        _cell_nodes.values[base + 2] = se
        _cell_nodes.values[base + 3] = ne
        _cell_nodes.values[base + 4] = nw
    end

    # cell → edges (4 per cell)
    _cell_edges = CSRMapping(n_cells, 4)
    for ilat in 1:nlat, ilon in 1:nlon
        cid = (ilat - 1) * nlon + ilon
        south = (ilat - 1) * nlon + ilon
        north = ilat * nlon + ilon
        west = n_h_edges + (ilat - 1) * (nlon + 1) + ilon
        east_ilon = ilon == nlon ? 1 : ilon + 1
        east = n_h_edges + (ilat - 1) * (nlon + 1) + east_ilon
        base = _cell_edges.offsets[cid] - 1
        _cell_edges.values[base + 1] = south
        _cell_edges.values[base + 2] = north
        _cell_edges.values[base + 3] = west
        _cell_edges.values[base + 4] = east
    end

    # edge → nodes (2 per edge)
    _edge_nodes = CSRMapping(n_edges, 2)
    # Horizontal edges
    for ilat in 1:(nlat + 1), ilon in 1:nlon
        eid = (ilat - 1) * nlon + ilon
        n1 = (ilat - 1) * (nlon + 1) + ilon
        n2 = (ilat - 1) * (nlon + 1) + ilon + 1
        base = _edge_nodes.offsets[eid] - 1
        _edge_nodes.values[base + 1] = n1
        _edge_nodes.values[base + 2] = n2
    end
    # Vertical edges
    for ilat in 1:nlat, ilon in 1:(nlon + 1)
        eid = n_h_edges + (ilat - 1) * (nlon + 1) + ilon
        n1 = (ilat - 1) * (nlon + 1) + ilon
        n2 = ilat * (nlon + 1) + ilon
        base = _edge_nodes.offsets[eid] - 1
        _edge_nodes.values[base + 1] = n1
        _edge_nodes.values[base + 2] = n2
    end

    # cell → cells (4 per cell, 0 sentinel at poles)
    _cell_cells = CSRMapping(n_cells, 4)
    for ilat in 1:nlat, ilon in 1:nlon
        cid = (ilat - 1) * nlon + ilon
        south = ilat > 1 ? (ilat - 2) * nlon + ilon : 0
        north = ilat < nlat ? ilat * nlon + ilon : 0
        west = ilon > 1 ? (ilat - 1) * nlon + (ilon - 1) : (ilat - 1) * nlon + nlon
        east = ilon < nlon ? (ilat - 1) * nlon + (ilon + 1) : (ilat - 1) * nlon + 1
        base = _cell_cells.offsets[cid] - 1
        _cell_cells.values[base + 1] = south
        _cell_cells.values[base + 2] = north
        _cell_cells.values[base + 3] = west
        _cell_cells.values[base + 4] = east
    end

    # edge → cells (variable: 1 at polar boundaries, 2 elsewhere)
    edge_cell_counts = fill(2, n_edges)
    # South pole horizontal edges: only 1 cell
    for ilon in 1:nlon
        eid = ilon
        edge_cell_counts[eid] = 1
    end
    # North pole horizontal edges: only 1 cell
    for ilon in 1:nlon
        eid = nlat * nlon + ilon
        edge_cell_counts[eid] = 1
    end
    _edge_cells, ptrs = CSRMapping(n_edges, edge_cell_counts)

    # Fill edge → cells
    for ilat in 1:nlat, ilon in 1:nlon
        cid = (ilat - 1) * nlon + ilon
        # south edge
        eid_s = (ilat - 1) * nlon + ilon
        _edge_cells.values[ptrs[eid_s]] = cid; ptrs[eid_s] += 1
        # north edge
        eid_n = ilat * nlon + ilon
        _edge_cells.values[ptrs[eid_n]] = cid; ptrs[eid_n] += 1
        # west edge
        eid_w = n_h_edges + (ilat - 1) * (nlon + 1) + ilon
        _edge_cells.values[ptrs[eid_w]] = cid; ptrs[eid_w] += 1
        # east edge
        east_ilon = ilon == nlon ? 1 : ilon + 1
        eid_e = n_h_edges + (ilat - 1) * (nlon + 1) + east_ilon
        _edge_cells.values[ptrs[eid_e]] = cid; ptrs[eid_e] += 1
    end

    # node → edges (variable: 2 at poles, 4 interior, 3 at boundaries)
    node_edge_counts = fill(4, n_nodes)
    # Corner nodes (south pole corners)
    for ilon in [1, nlon + 1]
        nid = ilon
        node_edge_counts[nid] = 2
    end
    # North pole corners
    for ilon in [1, nlon + 1]
        nid = nlat * (nlon + 1) + ilon
        node_edge_counts[nid] = 2
    end
    _node_edges, ptrs = CSRMapping(n_nodes, node_edge_counts)

    # Fill node → edges
    # Horizontal edges
    for ilat in 1:(nlat + 1), ilon in 1:nlon
        eid = (ilat - 1) * nlon + ilon
        n1 = (ilat - 1) * (nlon + 1) + ilon
        n2 = (ilat - 1) * (nlon + 1) + ilon + 1
        _node_edges.values[ptrs[n1]] = eid; ptrs[n1] += 1
        _node_edges.values[ptrs[n2]] = eid; ptrs[n2] += 1
    end
    # Vertical edges
    for ilat in 1:nlat, ilon in 1:(nlon + 1)
        eid = n_h_edges + (ilat - 1) * (nlon + 1) + ilon
        n1 = (ilat - 1) * (nlon + 1) + ilon
        n2 = ilat * (nlon + 1) + ilon
        _node_edges.values[ptrs[n1]] = eid; ptrs[n1] += 1
        _node_edges.values[ptrs[n2]] = eid; ptrs[n2] += 1
    end

    return LatLonGrid(
        M, lat_edges, lon_edges, R, nlat, nlon, nodes, cell_volumes, cell_centroids,
        _cell_nodes, _cell_edges, _cell_cells, _edge_nodes, _edge_cells, _node_edges,
        Ref{Union{Nothing, AbstractManifoldMesh{typeof(M)}}}(nothing))
end

function _cell_indices(g::LatLonGrid, cell_id::Int)
    ilat = div(cell_id - 1, g.nlon) + 1
    ilon = rem(cell_id - 1, g.nlon) + 1
    return (ilat, ilon)
end

function _cell_linear_index(g::LatLonGrid, ilat::Int, ilon::Int)
    return (ilat - 1) * g.nlon + ilon
end

function _node_indices(g::LatLonGrid, node_id::Int)
    ilat = div(node_id - 1, g.nlon + 1) + 1
    ilon = rem(node_id - 1, g.nlon + 1) + 1
    return (ilat, ilon)
end

function _node_linear_index(g::LatLonGrid, ilat::Int, ilon::Int)
    return (ilat - 1) * (g.nlon + 1) + ilon
end

# -- TopologyStyle Override --

TopologyStyle(::Type{<:LatLonGrid}) = IsGrid()
CellTypeStyle(::Type{<:LatLonGrid}) = IsUniform{4}()
PatchStyle(::Type{<:LatLonGrid}) = NoPatch()

has_dual(g::LatLonGrid) = g._dual[] !== nothing

# -- Global Information --

manifold(g::LatLonGrid) = g.manifold
num_cells(g::LatLonGrid) = g.nlat * g.nlon
num_nodes(g::LatLonGrid) = (g.nlat + 1) * (g.nlon + 1)
function num_edges(g::LatLonGrid)
    return (g.nlat + 1) * g.nlon + g.nlat * (g.nlon + 1)
end

# -- Geometry: node_coordinates --

function node_coordinates(g::LatLonGrid, node_id::Int)
    _check_node_id(g, node_id)
    ilat, ilon = _node_indices(g, node_id)
    return g.nodes[ilat, ilon]
end

# -- Geometry: cell_volume (cache read) --

function cell_volume(g::LatLonGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    ilat, ilon = _cell_indices(g, cell_id)
    return g.cell_volumes[ilat, ilon]
end

all_cell_volumes(g::LatLonGrid) = vec(g.cell_volumes)
all_node_coordinates(g::LatLonGrid) = vec(g.nodes)

# -- Geometry: cell_centroid (cache read) --

function cell_centroid(g::LatLonGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    ilat, ilon = _cell_indices(g, cell_id)
    return g.cell_centroids[ilat, ilon]
end

# -- Connectivity --

function cell_nodes(g::LatLonGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return getindex_fixed(g._cell_nodes, cell_id, Val(4))
end

function cell_cells(g::LatLonGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return getindex_fixed(g._cell_cells, cell_id, Val(4))
end

function node_cells(g::LatLonGrid, node_id::Int)
    _check_node_id(g, node_id)
    ilat, ilon = _node_indices(g, node_id)
    cells = Int[]
    for i in (ilat - 1, ilat)
        (i < 1 || i > g.nlat) && continue
        for j in (ilon - 1, ilon)
            jj = j < 1 ? g.nlon : (j > g.nlon ? 1 : j)
            if 1 ≤ i ≤ g.nlat && 1 ≤ jj ≤ g.nlon
                push!(cells, _cell_linear_index(g, i, jj))
            end
        end
    end
    return cells
end

function cell_edges(g::LatLonGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return getindex_fixed(g._cell_edges, cell_id, Val(4))
end

function edge_nodes(g::LatLonGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    return getindex_fixed(g._edge_nodes, edge_id, Val(2))
end

function edge_cells(g::LatLonGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    return g._edge_cells[edge_id]
end

function node_edges(g::LatLonGrid, node_id::Int)
    _check_node_id(g, node_id)
    return g._node_edges[node_id]
end

# -- Internal: Edge Endpoints --

function _edge_endpoints(g::LatLonGrid, edge_id::Int)
    n1, n2 = getindex_fixed(g._edge_nodes, edge_id, Val(2))
    return (node_coordinates(g, n1), node_coordinates(g, n2))
end

# -- Geometry: edge_length --

function edge_length(g::LatLonGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    n1, n2 = _edge_endpoints(g, edge_id)
    return Manifolds.distance(g.manifold, n1, n2)
end

# -- Geometry: edge_midpoint --

function edge_midpoint(g::LatLonGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    n1, n2 = _edge_endpoints(g, edge_id)
    if Manifolds.distance(g.manifold, n1, n2) < 1e-14
        return SVector{3, Float64}(n1)
    end
    return SVector{3, Float64}(Manifolds.mid_point(g.manifold, n1, n2))
end

# -- Geometry: edge_outward_normal --

function edge_outward_normal(g::LatLonGrid, edge_id::Int, cell_id::Int)
    _check_edge_id(g, edge_id)
    _check_cell_id(g, cell_id)
    n1, n2 = _edge_endpoints(g, edge_id)

    # Guard: polar collapse -- degenerate edge with coincident endpoints
    if Manifolds.distance(g.manifold, n1, n2) < 1e-14
        midpoint = n1
        return (base_point = SVector{3, Float64}(midpoint),
            normal = zero(SVector{3, Float64}))
    end

    midpoint = Manifolds.mid_point(g.manifold, n1, n2)
    gc_normal = cross(SVector(n1), SVector(n2))  # great circle plane normal

    tangent = normalize(cross(gc_normal, midpoint))
    cell_c = cell_centroid(g, cell_id)
    cell_side = sign(dot(gc_normal, SVector(cell_c)))

    outward = cell_side * cross(tangent, SVector(midpoint))
    outward = Manifolds.project(g.manifold, midpoint, outward)

    return (base_point = SVector{3, Float64}(midpoint),
        normal = SVector{3, Float64}(outward))
end

# -- Boundary Markers --

boundary_nodes(g::LatLonGrid, marker) = Int[]
boundary_edges(g::LatLonGrid, marker) = Int[]
