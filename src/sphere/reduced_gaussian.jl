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

struct ReducedGaussianGrid{M <: AbstractManifold} <: AbstractManifoldMesh{M}
    manifold::M
    nlat::Int
    R::Float64
    lat_points::Vector{Float64}     # Gaussian latitudes in radians
    lon_counts::Vector{Int}         # Number of longitude cells per latitude band
    nodes::Vector{SVector{3, Float64}}
    cell_volumes::Vector{Float64}
    cell_centroids::Vector{SVector{3, Float64}}
    _cell_nodes::Vector{NTuple{4, Int}}
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

# -- Constructor --

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

    return ReducedGaussianGrid{typeof(M)}(
        M, nlat, R, lat_points, lon_counts, nodes,
        cell_volumes, cell_centroids, _cell_nodes,
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

# -- Geometry --

function node_coordinates(g::ReducedGaussianGrid, node_id::Int)
    _check_node_id(g, node_id)
    return g.nodes[node_id]
end

function cell_volume(g::ReducedGaussianGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g.cell_volumes[cell_id]
end

function cell_centroid(g::ReducedGaussianGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g.cell_centroids[cell_id]
end

# -- Topology Stubs --

function cell_nodes(g::ReducedGaussianGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g._cell_nodes[cell_id]
end

# TODO(phase3): implement full topology (cell_cells, node_cells, cell_edges)
function cell_cells(g::ReducedGaussianGrid, cell_id::Int)
    error("not yet implemented")
end

function node_cells(g::ReducedGaussianGrid, node_id::Int)
    error("not yet implemented")
end

function cell_edges(g::ReducedGaussianGrid, cell_id::Int)
    error("not yet implemented")
end

# -- Edge Stubs --
# TODO(phase3): implement edge geometry

function edge_length(g::ReducedGaussianGrid, edge_id::Int)
    error("not yet implemented")
end

function edge_midpoint(g::ReducedGaussianGrid, edge_id::Int)
    error("not yet implemented")
end

function edge_outward_normal(g::ReducedGaussianGrid, edge_id::Int, cell_id::Int)
    error("not yet implemented")
end

# -- Boundary --

boundary_nodes(g::ReducedGaussianGrid, marker) = Int[]
boundary_edges(g::ReducedGaussianGrid, marker) = Int[]
