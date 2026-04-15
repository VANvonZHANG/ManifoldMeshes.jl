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
    cell_volumes = Matrix{Float64}(undef, nlat, nlon)
    for ilat in 1:nlat
        for ilon in 1:nlon
            ilon_next = ilon == nlon ? 1 : ilon + 1
            A = nodes[ilat, ilon]             # SW
            B = nodes[ilat, ilon_next]        # SE
            C = nodes[ilat + 1, ilon_next]    # NE
            D = nodes[ilat + 1, ilon]         # NW
            cell_volumes[ilat, ilon] =
                _spherical_triangle_area(R, A, B, C) +
                _spherical_triangle_area(R, A, C, D)
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
