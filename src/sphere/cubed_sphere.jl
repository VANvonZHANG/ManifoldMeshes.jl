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
        rotation = SMatrix{3, 3, Float64, 9}(I), R::Float64 = 1.0)
    n >= 1 || throw(ArgumentError("n must be >= 1, got $n"))
    projection in (:gnomonic, :equiangular) ||
        throw(ArgumentError("projection must be :gnomonic or :equiangular, got $projection"))
    R > 0 || throw(ArgumentError("R must be positive, got $R"))
    rotation = convert(SMatrix{3, 3, Float64, 9}, rotation)

    M = Sphere(2)

    # Each face has (n+1) x (n+1) nodes in local (s,t) coordinates.
    # For this initial implementation, nodes are NOT deduplicated across face
    # boundaries — each face stores its own copy of edge/corner nodes.
    nn_face = (n + 1) * (n + 1)
    face_node_offset = [0, 1, 2, 3, 4, 5] .* nn_face

    nodes = SVector{3, Float64}[]
    sizehint!(nodes, 6 * nn_face)

    for face in 1:6
        for j in 0:n
            t = -1.0 + 2.0 * j / n
            for i in 0:n
                s = -1.0 + 2.0 * i / n
                p = _cubed_sphere_face_point(face, s, t, projection)
                p = rotation * p
                p = R * normalize(p)
                push!(nodes, p)
            end
        end
    end

    # Pre-compute cell volumes (two spherical triangles per quad)
    ncells = 6 * n * n
    cell_volumes = Vector{Float64}(undef, ncells)
    cell_centroids = Vector{SVector{3, Float64}}(undef, ncells)

    for face in 1:6
        offset = face_node_offset[face]
        for j in 1:n
            for i in 1:n
                cell_id = _cubed_sphere_cell_id(n, face, i, j)

                # Node indices within face: row-major (i + 1, j + 1) for (n+1) x (n+1)
                sw = offset + (j - 1) * (n + 1) + i
                se = offset + (j - 1) * (n + 1) + (i + 1)
                ne = offset + j * (n + 1) + (i + 1)
                nw = offset + j * (n + 1) + i

                A = nodes[sw]
                B = nodes[se]
                C = nodes[ne]
                D = nodes[nw]

                area = spherical_triangle_area(R, A, B, C) +
                       spherical_triangle_area(R, A, C, D)
                cell_volumes[cell_id] = area

                verts = [A, B, C, D]
                c = Manifolds.mean(M, verts)
                cell_centroids[cell_id] = SVector{3, Float64}(c)
            end
        end
    end

    return CubedSphereGrid(
        M, n, R, nodes, cell_volumes, cell_centroids,
        Ref{Union{Nothing, AbstractManifoldMesh{typeof(M)}}}(nothing),
        rotation)
end

# -- Face Parameterization --

function _cubed_sphere_face_point(face::Int, s::Float64, t::Float64, projection::Symbol)
    if projection == :gnomonic
        return _gnomonic_point(face, s, t)
    else
        return _equiangular_point(face, s, t)
    end
end

function _gnomonic_point(face::Int, s::Float64, t::Float64)
    # Face normals:
    #   1: +Z, 2: -Z, 3: +Y, 4: -Y, 5: +X, 6: -X
    if face == 1
        p = SVector(s, t, 1.0)
    elseif face == 2
        p = SVector(-s, t, -1.0)
    elseif face == 3
        p = SVector(s, 1.0, -t)
    elseif face == 4
        p = SVector(s, -1.0, t)
    elseif face == 5
        p = SVector(1.0, t, -s)
    else  # face == 6
        p = SVector(-1.0, t, s)
    end
    return normalize(p)
end

function _equiangular_point(face::Int, s::Float64, t::Float64)
    a = atan(s)
    b = atan(t)
    cos_a = cos(a)
    sin_a = sin(a)
    cos_b = cos(b)
    sin_b = sin(b)

    if face == 1
        p = SVector(sin_a * cos_b, cos_a * sin_b, cos_a * cos_b)
    elseif face == 2
        p = SVector(-sin_a * cos_b, cos_a * sin_b, -cos_a * cos_b)
    elseif face == 3
        p = SVector(sin_a * cos_b, cos_a * cos_b, -cos_a * sin_b)
    elseif face == 4
        p = SVector(sin_a * cos_b, -cos_a * cos_b, cos_a * sin_b)
    elseif face == 5
        p = SVector(cos_a * cos_b, cos_a * sin_b, -sin_a * cos_b)
    else  # face == 6
        p = SVector(-cos_a * cos_b, cos_a * sin_b, sin_a * cos_b)
    end
    return normalize(p)
end

function _cubed_sphere_cell_id(n::Int, face::Int, i::Int, j::Int)
    return (face - 1) * n * n + (j - 1) * n + i
end

# -- Trait Implementations --

TopologyStyle(::Type{<:CubedSphereGrid}) = IsGrid()
CellTypeStyle(::Type{<:CubedSphereGrid}) = IsUniform{4}()
PatchStyle(::Type{<:CubedSphereGrid}) = MultiPatch{6}()

has_dual(g::CubedSphereGrid) = g._dual[] !== nothing

# -- Global Information --

manifold(g::CubedSphereGrid) = g.manifold
num_cells(g::CubedSphereGrid) = 6 * g.n * g.n

function num_nodes(g::CubedSphereGrid)
    n = g.n
    # Nodes are NOT deduplicated across face boundaries in this implementation.
    # Each face stores its own copy of edge/corner nodes.
    return 6 * (n + 1) * (n + 1)
end

function num_edges(g::CubedSphereGrid)
    n = g.n
    return 12 * n * (n + 1)
end

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
    face = div(cell_id - 1, n * n) + 1
    local_id = rem(cell_id - 1, n * n) + 1
    j = div(local_id - 1, n) + 1
    i = rem(local_id - 1, n) + 1

    nn_face = (n + 1) * (n + 1)
    offset = (face - 1) * nn_face

    sw = offset + (j - 1) * (n + 1) + i
    se = offset + (j - 1) * (n + 1) + (i + 1)
    ne = offset + j * (n + 1) + (i + 1)
    nw = offset + j * (n + 1) + i

    return (sw, se, ne, nw)
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
    n = g.n
    return div(cell_id - 1, n * n) + 1
end

function cell_local_2d(g::CubedSphereGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    n = g.n
    local_id = rem(cell_id - 1, n * n) + 1
    j = div(local_id - 1, n) + 1
    i = rem(local_id - 1, n) + 1
    return (i, j)
end
