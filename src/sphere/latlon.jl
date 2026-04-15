using Manifolds: Sphere

struct LatLonGrid{M<:AbstractManifold} <: AbstractManifoldMesh{M}
    manifold::M
    lat_edges::Vector{Float64}
    lon_edges::Vector{Float64}
    R::Float64
    nlat::Int
    nlon::Int
    nodes::Matrix{SVector{3,Float64}}
    cell_volumes::Matrix{Float64}
    cell_centroids::Matrix{SVector{3,Float64}}
end

function LatLonGrid(; lat_edges::Vector{Float64}, lon_edges::Vector{Float64}, R::Float64=1.0)
    # Validate
    length(lat_edges) < 2 && throw(ArgumentError("lat_edges must have at least 2 elements"))
    length(lon_edges) < 2 && throw(ArgumentError("lon_edges must have at least 2 elements"))
    lat_edges[1] != -90.0 && throw(ArgumentError("lat_edges must start at -90"))
    lat_edges[end] != 90.0 && throw(ArgumentError("lat_edges must end at 90"))
    !issorted(lat_edges) && throw(ArgumentError("lat_edges must be ascending"))
    lon_edges[1] != 0.0 && throw(ArgumentError("lon_edges must start at 0"))
    lon_edges[end] != 360.0 && throw(ArgumentError("lon_edges must end at 360"))
    !issorted(lon_edges) && throw(ArgumentError("lon_edges must be ascending"))

    nlat = length(lat_edges) - 1
    nlon = length(lon_edges) - 1
    M = Sphere(2)

    # Build nodes matrix [ilat, ilon], size (nlat+1) x (nlon+1)
    nodes = Matrix{SVector{3,Float64}}(undef, nlat + 1, nlon + 1)
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
            area_ac = _spherical_triangle_area(R, A, B, C) +
                       _spherical_triangle_area(R, A, C, D)
            if area_ac == 0.0
                # Degenerate A-C diagonal; use B-D diagonal instead
                area_ac = _spherical_triangle_area(R, B, C, D) +
                           _spherical_triangle_area(R, A, B, D)
            end
            if area_ac == 0.0
                # Both diagonals degenerate (e.g. 180° polar cells where
                # all vertices collapse to 2 antipodal points). Use lune formula.
                Δlon_rad = deg2rad(lon_edges[ilon_next] - lon_edges[ilon])
                area_ac = R^2 * (sin(deg2rad(lat_edges[ilat+1])) -
                                 sin(deg2rad(lat_edges[ilat]))) * abs(Δlon_rad)
            end
            cell_volumes[ilat, ilon] = area_ac
        end
    end

    # Pre-compute cell centroids via Manifolds.mean
    cell_centroids = Matrix{SVector{3,Float64}}(undef, nlat, nlon)
    for ilat in 1:nlat
        for ilon in 1:nlon
            ilon_next = ilon == nlon ? 1 : ilon + 1
            verts = [nodes[ilat, ilon], nodes[ilat, ilon_next],
                     nodes[ilat+1, ilon_next], nodes[ilat+1, ilon]]
            c = Manifolds.mean(M, verts)
            cell_centroids[ilat, ilon] = SVector{3,Float64}(c)
        end
    end

    return LatLonGrid(M, lat_edges, lon_edges, R, nlat, nlon, nodes, cell_volumes, cell_centroids)
end

# -- Internal: Spherical Triangle Area (l'Huilier's formula) --

function _spherical_triangle_area(R::Float64, A::SVector{3,Float64},
                                  B::SVector{3,Float64}, C::SVector{3,Float64})
    a = acos(clamp(dot(B, C) / (R * R), -1, 1))
    b = acos(clamp(dot(A, C) / (R * R), -1, 1))
    c = acos(clamp(dot(A, B) / (R * R), -1, 1))
    s = (a + b + c) / 2
    tan_half = tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2)
    E = 4 * atan(sqrt(max(tan_half, 0.0)))
    return R^2 * E
end

# -- Internal: Index Helpers --

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

# -- Global Information --

manifold(g::LatLonGrid) = g.manifold
num_cells(g::LatLonGrid) = g.nlat * g.nlon
num_nodes(g::LatLonGrid) = (g.nlat + 1) * (g.nlon + 1)
function num_edges(g::LatLonGrid)
    return (g.nlat + 1) * g.nlon + g.nlat * g.nlon
end

# -- Geometry: node_coordinates --

function node_coordinates(g::LatLonGrid, node_id::Int)
    ilat, ilon = _node_indices(g, node_id)
    return g.nodes[ilat, ilon]
end

# -- Geometry: cell_volume (cache read) --

function cell_volume(g::LatLonGrid, cell_id::Int)
    ilat, ilon = _cell_indices(g, cell_id)
    return g.cell_volumes[ilat, ilon]
end

# -- Geometry: cell_centroid (cache read) --

function cell_centroid(g::LatLonGrid, cell_id::Int)
    ilat, ilon = _cell_indices(g, cell_id)
    return g.cell_centroids[ilat, ilon]
end

# -- Connectivity --

function cell_nodes(g::LatLonGrid, cell_id::Int)
    ilat, ilon = _cell_indices(g, cell_id)
    ilon_next = ilon == g.nlon ? 1 : ilon + 1
    return (
        _node_linear_index(g, ilat,   ilon),       # SW
        _node_linear_index(g, ilat,   ilon_next),   # SE
        _node_linear_index(g, ilat+1, ilon_next),   # NE
        _node_linear_index(g, ilat+1, ilon),         # NW
    )
end

function cell_cells(g::LatLonGrid, cell_id::Int)
    ilat, ilon = _cell_indices(g, cell_id)
    south = ilat > 1     ? _cell_linear_index(g, ilat - 1, ilon) : 0
    north = ilat < g.nlat ? _cell_linear_index(g, ilat + 1, ilon) : 0
    west  = _cell_linear_index(g, ilat, ilon == 1    ? g.nlon : ilon - 1)
    east  = _cell_linear_index(g, ilat, ilon == g.nlon ? 1      : ilon + 1)
    return (south, north, west, east)
end

function node_cells(g::LatLonGrid, node_id::Int)
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
    ilat, ilon = _cell_indices(g, cell_id)
    n_h = (g.nlat + 1) * g.nlon
    south = (ilat - 1) * g.nlon + ilon
    north = ilat * g.nlon + ilon
    west  = n_h + (ilat - 1) * g.nlon + ilon
    east_ilon = ilon == g.nlon ? 1 : ilon + 1
    east  = n_h + (ilat - 1) * g.nlon + east_ilon
    return (south, north, west, east)
end

# -- Internal: Edge Endpoints --

function _edge_endpoints(g::LatLonGrid, edge_id::Int)
    n_h = (g.nlat + 1) * g.nlon
    if edge_id <= n_h
        # Horizontal edge: along a latitude circle
        idx = edge_id - 1
        ilat = div(idx, g.nlon) + 1
        ilon = rem(idx, g.nlon) + 1
        ilon_next = ilon == g.nlon ? 1 : ilon + 1
        return (g.nodes[ilat, ilon], g.nodes[ilat, ilon_next])
    else
        # Vertical edge: along a longitude line
        idx = edge_id - n_h - 1
        ilat = div(idx, g.nlon) + 1
        ilon = rem(idx, g.nlon) + 1
        return (g.nodes[ilat, ilon], g.nodes[ilat + 1, ilon])
    end
end

# -- Geometry: edge_length --

function edge_length(g::LatLonGrid, edge_id::Int)
    n1, n2 = _edge_endpoints(g, edge_id)
    return Manifolds.distance(g.manifold, n1, n2)
end

# -- Geometry: edge_midpoint --

function edge_midpoint(g::LatLonGrid, edge_id::Int)
    n1, n2 = _edge_endpoints(g, edge_id)
    if Manifolds.distance(g.manifold, n1, n2) < 1e-14
        return SVector{3,Float64}(n1)
    end
    return SVector{3,Float64}(Manifolds.mid_point(g.manifold, n1, n2))
end

# -- Geometry: edge_outward_normal --

function edge_outward_normal(g::LatLonGrid, edge_id::Int, cell_id::Int)
    n1, n2 = _edge_endpoints(g, edge_id)

    # Guard: polar collapse -- degenerate edge with coincident endpoints
    if Manifolds.distance(g.manifold, n1, n2) < 1e-14
        midpoint = n1
        return (base_point = SVector{3,Float64}(midpoint),
                normal = zero(SVector{3,Float64}))
    end

    midpoint = Manifolds.mid_point(g.manifold, n1, n2)
    gc_normal = cross(SVector(n1), SVector(n2))  # great circle plane normal

    tangent = normalize(cross(gc_normal, midpoint))
    cell_c = cell_centroid(g, cell_id)
    cell_side = sign(dot(gc_normal, SVector(cell_c)))

    outward = cell_side * cross(tangent, SVector(midpoint))
    outward = Manifolds.project(g.manifold, midpoint, outward)

    return (base_point = SVector{3,Float64}(midpoint),
            normal = SVector{3,Float64}(outward))
end

# -- Boundary Markers --

boundary_nodes(g::LatLonGrid, marker) = Int[]
boundary_edges(g::LatLonGrid, marker) = Int[]
