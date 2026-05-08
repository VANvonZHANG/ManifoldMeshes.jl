"""
    slerp(p1::SVector{3,Float64}, p2::SVector{3,Float64}, n::Int) -> Vector{Point3f}

Spherical linear interpolation between `p1` and `p2`, producing `n` evenly spaced
points on the great circle arc. For coincident endpoints (θ < 1e-10), returns a
single-point vector. For antipodal endpoints (θ ≈ π), falls back to a stable
intermediate great circle plane.
"""
function slerp(p1::SVector{3, Float64}, p2::SVector{3, Float64}, n::Int)
    d = clamp(dot(p1, p2), -1.0, 1.0)
    θ = acos(d)

    # Coincident endpoints
    if θ < 1e-10
        return [Point3f(Float32.(p1))]
    end

    sinθ = sin(θ)

    # Antipodal endpoints: pick an arbitrary perpendicular axis for the great circle
    if θ > π - 1e-10
        # Find a vector not parallel to p1
        ref = abs(p1[1]) < 0.9 ? SVector(1.0, 0.0, 0.0) : SVector(0.0, 1.0, 0.0)
        perp = normalize(cross(p1, ref))
        return [Point3f(Float32.(cos(t * π) * p1 + sin(t * π) * perp))
                for t in LinRange(0, 1, n)]
    end

    return [Point3f(Float32.(
        sin((1 - t) * θ) / sinθ * p1[i] + sin(t * θ) / sinθ * p2[i]
    )) for t in LinRange(0, 1, n) for i in 1:3]
end

"""
    node_points(g::AbstractManifoldMesh) -> Vector{Point3f}

Return all node positions as `Point3f` values, indexed by node ID.
"""
function node_points(g::AbstractManifoldMesh)
    return [Point3f(Float32.(node_coordinates(g, i))) for i in 1:num_nodes(g)]
end

"""
    _build_edge_endpoint_map(g::AbstractManifoldMesh) -> Dict{Int, Tuple{Int,Int}}

Reconstruct edge-to-endpoint-node mapping from the public cell interface.
For each cell, edges are ordered (south, north, west, east) and nodes are
ordered (SW, SE, NE, NW), giving:
  - south edge: SW -> SE  (nodes 1-2)
  - north edge: NW -> NE  (nodes 4-3)
  - west edge:  SW -> NW  (nodes 1-4)
  - east edge:  SE -> NE  (nodes 2-3)
"""
function _build_edge_endpoint_map(g::AbstractManifoldMesh)
    # node_pair_for_edge[edge_id] = (node_a, node_b)
    node_pair_for_edge = Dict{Int, Tuple{Int, Int}}()
    for cid in 1:num_cells(g)
        ns = cell_nodes(g, cid)
        es = cell_edges(g, cid)
        # south edge: SW(1) -> SE(2)
        _register_edge_pair!(node_pair_for_edge, es[1], ns[1], ns[2])
        # north edge: NW(4) -> NE(3)
        _register_edge_pair!(node_pair_for_edge, es[2], ns[4], ns[3])
        # west edge: SW(1) -> NW(4)
        _register_edge_pair!(node_pair_for_edge, es[3], ns[1], ns[4])
        # east edge: SE(2) -> NE(3)
        _register_edge_pair!(node_pair_for_edge, es[4], ns[2], ns[3])
    end
    return node_pair_for_edge
end

function _register_edge_pair!(
        d::Dict{Int, Tuple{Int, Int}}, edge_id::Int, a::Int, b::Int)
    if !haskey(d, edge_id)
        d[edge_id] = (a, b)
    end
end

"""
    edge_segments(g::AbstractManifoldMesh; n_arc_points::Int=20) -> Vector{Vector{Point3f}}

Return discretized great-circle arcs for each edge, as a vector of `Point3f` arrays.
Degenerate edges (coincident endpoints) return a single-point segment.
"""
function edge_segments(g::AbstractManifoldMesh; n_arc_points::Int=20)
    edge_map = _build_edge_endpoint_map(g)
    segments = Vector{Vector{Point3f}}(undef, num_edges(g))
    for eid in 1:num_edges(g)
        n1_id, n2_id = edge_map[eid]
        p1 = node_coordinates(g, n1_id)
        p2 = node_coordinates(g, n2_id)
        segments[eid] = slerp(p1, p2, n_arc_points)
    end
    return segments
end

"""
    cell_polygons(g::AbstractManifoldMesh; n_arc_points::Int=20) -> Vector{Vector{Point3f}}

Return closed cell boundaries as discretized great-circle polygons.
Each polygon is a vector of `Point3f` where the first point equals the last.
Edges are ordered: south, east, north (reversed), west (reversed).
"""
function cell_polygons(g::AbstractManifoldMesh; n_arc_points::Int=20)
    edge_map = _build_edge_endpoint_map(g)
    polygons = Vector{Vector{Point3f}}(undef, num_cells(g))

    for cid in 1:num_cells(g)
        ns = cell_nodes(g, cid)
        es = cell_edges(g, cid)

        # Build arcs for each edge
        south_arc = slerp(node_coordinates(g, ns[1]), node_coordinates(g, ns[2]), n_arc_points)
        east_arc = slerp(node_coordinates(g, ns[2]), node_coordinates(g, ns[3]), n_arc_points)
        north_arc = slerp(node_coordinates(g, ns[4]), node_coordinates(g, ns[3]), n_arc_points)
        west_arc = slerp(node_coordinates(g, ns[1]), node_coordinates(g, ns[4]), n_arc_points)

        # Concatenate: south + east + north(reversed) + west(reversed)
        # Skip first point of each subsequent arc to avoid duplication at corners
        total_len = length(south_arc) + length(east_arc) - 1 +
                    length(north_arc) - 1 + length(west_arc) - 1
        polygon = Vector{Point3f}(undef, total_len)

        idx = 1
        for p in south_arc
            polygon[idx] = p; idx += 1
        end
        for i in 2:length(east_arc)
            polygon[idx] = east_arc[i]; idx += 1
        end
        for i in 2:length(north_arc)
            polygon[idx] = north_arc[length(north_arc) - i + 2]; idx += 1
        end
        for i in 2:length(west_arc)
            polygon[idx] = west_arc[length(west_arc) - i + 2]; idx += 1
        end

        polygons[cid] = polygon
    end

    return polygons
end
