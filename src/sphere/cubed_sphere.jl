using Manifolds: Sphere
using StaticArrays: SMatrix

struct CubedSphereGrid{M <: AbstractManifold} <: AbstractManifoldMesh{M}
    manifold::M
    n::Int
    R::Float64
    nodes::Vector{SVector{3, Float64}}
    cell_volumes::Vector{Float64}
    cell_centroids::Vector{SVector{3, Float64}}
    _dual::Base.RefValue{Union{Nothing, AbstractManifoldMesh{M}}}
    rotation::SMatrix{3, 3, Float64, 9}
    face_local_nodes::Vector{Matrix{Int}}  # maps (face, j, i) -> global node id
end

# -- Internal: Bounds Checking --

@inline function _check_cell_id(g::CubedSphereGrid, cell_id::Int)
    @boundscheck 1 <= cell_id <= num_cells(g) ||
                 throw(BoundsError("cell_id $cell_id out of range [1, $(num_cells(g))]"))
    nothing
end

@inline function _check_node_id(g::CubedSphereGrid, node_id::Int)
    @boundscheck 1 <= node_id <= num_nodes(g) ||
                 throw(BoundsError("node_id $node_id out of range [1, $(num_nodes(g))]"))
    nothing
end

@inline function _check_edge_id(g::CubedSphereGrid, edge_id::Int)
    @boundscheck 1 <= edge_id <= num_edges(g) ||
                 throw(BoundsError("edge_id $edge_id out of range [1, $(num_edges(g))]"))
    nothing
end

function CubedSphereGrid(; n::Int, projection::Symbol = :gnomonic,
        rotation::SMatrix{3, 3, Float64, 9} = SMatrix{3, 3, Float64}(I),
        R::Float64 = 1.0)
    # Validate
    n < 1 && throw(ArgumentError("n must be >= 1, got $n"))
    projection ∉ (:gnomonic, :equiangular) &&
        throw(ArgumentError("projection must be :gnomonic or :equiangular, got $projection"))
    R <= 0 && throw(ArgumentError("R must be positive, got $R"))

    M = Sphere(2)

    # --- Generate face-local nodes and map to global unique nodes ---
    # Use a dictionary keyed by rounded coordinates to find duplicates
    tol = 1e-12
    node_map = Dict{SVector{3, Float64}, Int}()
    nodes = SVector{3, Float64}[]

    # face_local_nodes[face][j,i] = global_node_id
    face_local_nodes = [Matrix{Int}(undef, n + 1, n + 1) for _ in 1:6]

    for face in 1:6
        for j in 0:n, i in 0:n
            s = -1.0 + 2.0 * i / n
            t = -1.0 + 2.0 * j / n
            p = _cubed_sphere_face_point(face, s, t, projection)
            p = rotation * p
            p = R * normalize(p)

            # Round to tol to catch floating-point duplicates
            key = SVector(
                round(p[1] / tol) * tol,
                round(p[2] / tol) * tol,
                round(p[3] / tol) * tol)

            if haskey(node_map, key)
                face_local_nodes[face][j + 1, i + 1] = node_map[key]
            else
                global_id = length(nodes) + 1
                push!(nodes, SVector{3, Float64}(p))
                node_map[key] = global_id
                face_local_nodes[face][j + 1, i + 1] = global_id
            end
        end
    end

    # --- Pre-compute cell volumes and centroids ---
    n_cells = 6 * n * n
    cell_volumes = Vector{Float64}(undef, n_cells)
    cell_centroids = Vector{SVector{3, Float64}}(undef, n_cells)

    for face in 1:6
        for j in 1:n, i in 1:n
            cell_id = _cubed_sphere_cell_id(n, face, i, j)
            local_n = face_local_nodes[face]

            A = nodes[local_n[j, i]]       # SW
            B = nodes[local_n[j, i + 1]]   # SE
            C = nodes[local_n[j + 1, i + 1]]  # NE
            D = nodes[local_n[j + 1, i]]      # NW

            area = spherical_triangle_area(R, A, B, C) +
                   spherical_triangle_area(R, A, C, D)
            cell_volumes[cell_id] = area

            c = Manifolds.mean(M, (A, B, C, D))
            cell_centroids[cell_id] = SVector{3, Float64}(c)
        end
    end

    return CubedSphereGrid{typeof(M)}(
        M, n, R, nodes, cell_volumes, cell_centroids,
        Ref{Union{Nothing, AbstractManifoldMesh{typeof(M)}}}(nothing),
        rotation, face_local_nodes)
end

# -- Face Parameterization Functions --

function _cubed_sphere_face_point(face::Int, s::Float64, t::Float64, projection::Symbol)
    if projection === :gnomonic
        return _gnomonic_point(face, s, t)
    else
        return _equiangular_point(face, s, t)
    end
end

function _gnomonic_point(face::Int, s::Float64, t::Float64)
    # Face normals: +Z(1), -Z(2), +Y(3), -Y(4), +X(5), -X(6)
    if face == 1
        return normalize(SVector(s, t, 1.0))
    elseif face == 2
        return normalize(SVector(s, t, -1.0))
    elseif face == 3
        return normalize(SVector(s, 1.0, t))
    elseif face == 4
        return normalize(SVector(s, -1.0, t))
    elseif face == 5
        return normalize(SVector(1.0, s, t))
    else  # face == 6
        return normalize(SVector(-1.0, s, t))
    end
end

function _equiangular_point(face::Int, s::Float64, t::Float64)
    # Map s, t in [-1, 1] via atan to get angular coordinates
    α = atan(s)
    β = atan(t)
    if face == 1
        return normalize(SVector(tan(α), tan(β), 1.0))
    elseif face == 2
        return normalize(SVector(tan(α), tan(β), -1.0))
    elseif face == 3
        return normalize(SVector(tan(α), 1.0, tan(β)))
    elseif face == 4
        return normalize(SVector(tan(α), -1.0, tan(β)))
    elseif face == 5
        return normalize(SVector(1.0, tan(α), tan(β)))
    else  # face == 6
        return normalize(SVector(-1.0, tan(α), tan(β)))
    end
end

function _cubed_sphere_cell_id(n::Int, face::Int, i::Int, j::Int)
    return (face - 1) * n * n + (j - 1) * n + i
end

# -- Trait Implementations --

TopologyStyle(::Type{<:CubedSphereGrid}) = IsGrid()
CellTypeStyle(::Type{<:CubedSphereGrid}) = IsUniform{4}()
PatchStyle(::Type{<:CubedSphereGrid}) = MultiPatch{6}()

has_dual(g::CubedSphereGrid) = g._dual[] !== nothing

function dual(g::CubedSphereGrid)
    error("dual construction for CubedSphereGrid is not yet implemented")
end

# -- Global Information --

manifold(g::CubedSphereGrid) = g.manifold
num_cells(g::CubedSphereGrid) = 6 * g.n * g.n
num_nodes(g::CubedSphereGrid) = length(g.nodes)
num_edges(g::CubedSphereGrid) = 12 * g.n * (g.n + 1)

# -- Geometry (with @boundscheck) --

function node_coordinates(g::CubedSphereGrid, node_id::Int)
    _check_node_id(g, node_id)
    return g.nodes[node_id]
end

function cell_volume(g::CubedSphereGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g.cell_volumes[cell_id]
end

function cell_centroid(g::CubedSphereGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g.cell_centroids[cell_id]
end

# -- Topology --

function cell_nodes(g::CubedSphereGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    n = g.n
    face, (i, j) = cell_face(g, cell_id), cell_local_2d(g, cell_id)
    fln = g.face_local_nodes[face]
    return (
        fln[j, i],       # SW
        fln[j, i + 1],   # SE
        fln[j + 1, i + 1],  # NE
        fln[j + 1, i]       # NW
    )
end

function cell_cells(g::CubedSphereGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    error("cell_cells not yet implemented for CubedSphereGrid")
end

function node_cells(g::CubedSphereGrid, node_id::Int)
    _check_node_id(g, node_id)
    error("node_cells not yet implemented for CubedSphereGrid")
end

function cell_edges(g::CubedSphereGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    error("cell_edges not yet implemented for CubedSphereGrid")
end

# -- Edge Stubs --

function edge_length(g::CubedSphereGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    error("edge_length not yet implemented for CubedSphereGrid")
end

function edge_midpoint(g::CubedSphereGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    error("edge_midpoint not yet implemented for CubedSphereGrid")
end

function edge_outward_normal(g::CubedSphereGrid, edge_id::Int, cell_id::Int)
    _check_edge_id(g, edge_id)
    _check_cell_id(g, cell_id)
    error("edge_outward_normal not yet implemented for CubedSphereGrid")
end

# -- Boundary --

boundary_nodes(g::CubedSphereGrid, marker) = Int[]
boundary_edges(g::CubedSphereGrid, marker) = Int[]

# -- Patch Queries --

function cell_face(g::CubedSphereGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return div(cell_id - 1, g.n * g.n) + 1
end

function cell_local_2d(g::CubedSphereGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    n = g.n
    local_id = rem(cell_id - 1, n * n) + 1
    i = rem(local_id - 1, n) + 1
    j = div(local_id - 1, n) + 1
    return (i, j)
end
