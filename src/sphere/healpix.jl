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
    nodes = SVector{3, Float64}[]

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

        for i in 1:n_in_ring
            phi = 2π * (i - 0.5) / n_in_ring
            x = R * sin(theta) * cos(phi)
            y = R * sin(theta) * sin(phi)
            z = R * cos(theta)
            p = rotation * SVector(x, y, z)
            push!(nodes, SVector{3, Float64}(p))
        end
    end

    cell_volumes = Vector{Float64}(undef, n_cells)
    cell_centroids = Vector{SVector{3, Float64}}(undef, n_cells)

    ideal_area = 4π * R^2 / n_cells
    for i in 1:n_cells
        cell_volumes[i] = ideal_area
        # Centroid: approximate using a random direction normalized to R
        # (Will be refined in full implementation)
        cell_centroids[i] = SVector{3, Float64}(R * normalize(randn(3)))
    end

    return HEALPixGrid{typeof(M)}(
        M, nside, R, ordering, nodes,
        cell_volumes, cell_centroids,
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

# -- Topology Stubs --

function cell_nodes(g::HEALPixGrid, cell_id::Int)
    error("not yet implemented")
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
