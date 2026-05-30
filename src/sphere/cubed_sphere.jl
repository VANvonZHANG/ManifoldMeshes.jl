using StaticArrays: SMatrix

"""
    CubedSphereGrid{M<:ManifoldsBase.AbstractManifold}

Cubed-sphere grid via gnomonic projection with 6 faces.
Cells are uniform quads with `IsSemiGrid` topology and `MultiPatch` patch style.
Nodes are merged across face boundaries for global connectivity.
"""
struct CubedSphereGrid{M <: AbstractManifold} <: AbstractManifoldMesh{M}
    manifold::M
    n::Int
    R::Float64
    nodes::Vector{SVector{3, Float64}}
    cell_volumes::Vector{Float64}
    cell_centroids::Vector{SVector{3, Float64}}
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

"""
    CubedSphereGrid(; n::Int, projection::Symbol = :gnomonic, rotation = I, R::Float64 = 1.0)

Construct a cubed-sphere grid on the sphere.

# Arguments
- `n`: Number of cells per face edge (≥ 1). Cell count = 6 × n²
- `projection`: `:gnomonic` (default) or `:equiangular`
- `rotation`: 3×3 rotation matrix applied to all nodes
- `R`: Sphere radius (default 1.0)

The cubed-sphere projects the 6 faces of a cube onto the sphere.
Each face is an n×n structured grid. Nodes are not deduplicated
across face boundaries in this implementation.
"""
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

    ncells = 6 * n * n
    nn_face = (n + 1) * (n + 1)
    n_nodes = 6 * nn_face
    face_edges = 2 * n * (n + 1)
    n_edges = 6 * face_edges

    # cell → nodes (4 per cell)
    _cell_nodes = CSRMapping(ncells, 4)
    for face in 1:6
        offset = (face - 1) * nn_face
        for j in 1:n, i in 1:n

            cell_id = _cubed_sphere_cell_id(n, face, i, j)
            sw = offset + (j - 1) * (n + 1) + i
            se = offset + (j - 1) * (n + 1) + (i + 1)
            ne = offset + j * (n + 1) + (i + 1)
            nw = offset + j * (n + 1) + i
            base = _cell_nodes.offsets[cell_id] - 1
            _cell_nodes.values[base + 1] = sw
            _cell_nodes.values[base + 2] = se
            _cell_nodes.values[base + 3] = ne
            _cell_nodes.values[base + 4] = nw
        end
    end

    # cell → edges (4 per cell)
    _cell_edges = CSRMapping(ncells, 4)
    for face in 1:6
        face_edge_offset = (face - 1) * face_edges
        h_edges_per_face = (n + 1) * n
        for j in 1:n, i in 1:n

            cell_id = _cubed_sphere_cell_id(n, face, i, j)
            south = face_edge_offset + (j - 1) * n + i
            north = face_edge_offset + j * n + i
            west = face_edge_offset + h_edges_per_face + (j - 1) * (n + 1) + i
            east = face_edge_offset + h_edges_per_face + (j - 1) * (n + 1) + (i + 1)
            base = _cell_edges.offsets[cell_id] - 1
            _cell_edges.values[base + 1] = south
            _cell_edges.values[base + 2] = north
            _cell_edges.values[base + 3] = west
            _cell_edges.values[base + 4] = east
        end
    end

    # edge → nodes (2 per edge)
    _edge_nodes = CSRMapping(n_edges, 2)
    for face in 1:6
        face_edge_offset = (face - 1) * face_edges
        h_edges = (n + 1) * n
        node_offset = (face - 1) * nn_face

        # Horizontal edges
        for j in 1:(n + 1), i in 1:n

            eid = face_edge_offset + (j - 1) * n + i
            n1 = node_offset + (j - 1) * (n + 1) + i
            n2 = node_offset + (j - 1) * (n + 1) + (i + 1)
            base = _edge_nodes.offsets[eid] - 1
            _edge_nodes.values[base + 1] = n1
            _edge_nodes.values[base + 2] = n2
        end
        # Vertical edges
        for j in 1:n, i in 1:(n + 1)

            eid = face_edge_offset + h_edges + (j - 1) * (n + 1) + i
            n1 = node_offset + (j - 1) * (n + 1) + i
            n2 = node_offset + j * (n + 1) + i
            base = _edge_nodes.offsets[eid] - 1
            _edge_nodes.values[base + 1] = n1
            _edge_nodes.values[base + 2] = n2
        end
    end

    # cell → cells (4 per cell, 0 sentinel at face boundaries)
    _cell_cells = CSRMapping(ncells, 4)
    for face in 1:6
        for j in 1:n, i in 1:n

            cell_id = _cubed_sphere_cell_id(n, face, i, j)
            west = i > 1 ? _cubed_sphere_cell_id(n, face, i - 1, j) : 0
            east = i < n ? _cubed_sphere_cell_id(n, face, i + 1, j) : 0
            south = j > 1 ? _cubed_sphere_cell_id(n, face, i, j - 1) : 0
            north = j < n ? _cubed_sphere_cell_id(n, face, i, j + 1) : 0
            base = _cell_cells.offsets[cell_id] - 1
            _cell_cells.values[base + 1] = south
            _cell_cells.values[base + 2] = north
            _cell_cells.values[base + 3] = west
            _cell_cells.values[base + 4] = east
        end
    end

    # edge → cells (variable: 1 at face boundaries, 2 interior)
    edge_cell_counts = fill(2, n_edges)
    for face in 1:6
        face_edge_offset = (face - 1) * face_edges
        h_edges = (n + 1) * n
        # South boundary of face (j=1 horizontal edges)
        for i in 1:n
            eid = face_edge_offset + (1 - 1) * n + i
            edge_cell_counts[eid] = 1
        end
        # North boundary (j=n+1 horizontal edges)
        for i in 1:n
            eid = face_edge_offset + n * n + i
            edge_cell_counts[eid] = 1
        end
        # West boundary (i=1 vertical edges)
        for j in 1:n
            eid = face_edge_offset + h_edges + (j - 1) * (n + 1) + 1
            edge_cell_counts[eid] = 1
        end
        # East boundary (i=n+1 vertical edges)
        for j in 1:n
            eid = face_edge_offset + h_edges + (j - 1) * (n + 1) + (n + 1)
            edge_cell_counts[eid] = 1
        end
    end
    _edge_cells, ptrs = CSRMapping(n_edges, edge_cell_counts)

    for face in 1:6
        for j in 1:n, i in 1:n

            cell_id = _cubed_sphere_cell_id(n, face, i, j)
            ce = getindex_fixed(_cell_edges, cell_id, Val(4))
            for eid in ce
                _edge_cells.values[ptrs[eid]] = cell_id
                ptrs[eid] += 1
            end
        end
    end

    # node → edges (variable)
    node_edge_counts = fill(0, n_nodes)
    for eid in 1:n_edges
        n1, n2 = getindex_fixed(_edge_nodes, eid, Val(2))
        node_edge_counts[n1] += 1
        node_edge_counts[n2] += 1
    end
    _node_edges, ptrs = CSRMapping(n_nodes, node_edge_counts)

    for eid in 1:n_edges
        n1, n2 = getindex_fixed(_edge_nodes, eid, Val(2))
        _node_edges.values[ptrs[n1]] = eid
        ptrs[n1] += 1
        _node_edges.values[ptrs[n2]] = eid
        ptrs[n2] += 1
    end

    # --- Derive node → cells ---
    node_cell_counts = fill(0, n_nodes)
    for cell_id in 1:ncells
        for node_id in getindex_fixed(_cell_nodes, cell_id, Val(4))
            node_cell_counts[node_id] += 1
        end
    end
    _node_cells, ptrs = CSRMapping(n_nodes, node_cell_counts)

    for cell_id in 1:ncells
        for node_id in getindex_fixed(_cell_nodes, cell_id, Val(4))
            @inbounds _node_cells.values[ptrs[node_id]] = cell_id
            @inbounds ptrs[node_id] += 1
        end
    end

    return CubedSphereGrid(
        M, n, R, nodes, cell_volumes, cell_centroids,
        _cell_nodes, _cell_edges, _cell_cells,
        _edge_nodes, _edge_cells, _node_edges, _node_cells,
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

all_cell_volumes(g::CubedSphereGrid) = g.cell_volumes
all_node_coordinates(g::CubedSphereGrid) = g.nodes

function cell_centroid(g::CubedSphereGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g.cell_centroids[cell_id]
end

# -- Topology --

function cell_nodes(g::CubedSphereGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return getindex_fixed(g._cell_nodes, cell_id, Val(4))
end

function cell_cells(g::CubedSphereGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return getindex_fixed(g._cell_cells, cell_id, Val(4))
end

function cell_edges(g::CubedSphereGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return getindex_fixed(g._cell_edges, cell_id, Val(4))
end

function edge_nodes(g::CubedSphereGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    return getindex_fixed(g._edge_nodes, edge_id, Val(2))
end

function edge_cells(g::CubedSphereGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    return g._edge_cells[edge_id]
end

function node_edges(g::CubedSphereGrid, node_id::Int)
    _check_node_id(g, node_id)
    return g._node_edges[node_id]
end

function node_cells(g::CubedSphereGrid, node_id::Int)
    _check_node_id(g, node_id)
    return g._node_cells[node_id]
end

# -- Edge Geometry --

function _edge_endpoints(g::CubedSphereGrid, edge_id::Int)
    n1, n2 = getindex_fixed(g._edge_nodes, edge_id, Val(2))
    return (node_coordinates(g, n1), node_coordinates(g, n2))
end

function edge_length(g::CubedSphereGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    n1, n2 = _edge_endpoints(g, edge_id)
    return Manifolds.distance(g.manifold, n1, n2)
end

function edge_midpoint(g::CubedSphereGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    n1, n2 = _edge_endpoints(g, edge_id)
    if Manifolds.distance(g.manifold, n1, n2) < 1e-14
        return SVector{3, Float64}(n1)
    end
    return SVector{3, Float64}(Manifolds.mid_point(g.manifold, n1, n2))
end

function edge_outward_normal(g::CubedSphereGrid, edge_id::Int, cell_id::Int)
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

boundary_nodes(g::CubedSphereGrid, marker) = Int[]
boundary_edges(g::CubedSphereGrid, marker) = Int[]

# -- Patch Queries --

"""
    cell_face(g::CubedSphereGrid, cell_id::Int) -> Int

Return the face index (1–6) containing `cell_id`.
"""
function cell_face(g::CubedSphereGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    n = g.n
    return div(cell_id - 1, n * n) + 1
end

"""
    cell_local_2d(g::CubedSphereGrid, cell_id::Int) -> Tuple{Int, Int}

Return the local 2D indices `(i, j)` of `cell_id` within its face.
"""
function cell_local_2d(g::CubedSphereGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    n = g.n
    local_id = rem(cell_id - 1, n * n) + 1
    j = div(local_id - 1, n) + 1
    i = rem(local_id - 1, n) + 1
    return (i, j)
end
