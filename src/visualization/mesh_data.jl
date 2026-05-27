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
                sin((1 - t) * θ) / sinθ * p1 .+ sin(t * θ) / sinθ * p2
            )) for t in LinRange(0, 1, n)]
end

"""
    node_points(g::AbstractManifoldMesh) -> Vector{Point3f}

Return all node positions as `Point3f` values, indexed by node ID.
"""
function node_points(g::AbstractManifoldMesh)
    return [Point3f(Float32.(node_coordinates(g, i))) for i in 1:num_nodes(g)]
end

"""
    edge_segments(g::AbstractManifoldMesh; n_arc_points::Int=20) -> Vector{Vector{Point3f}}

Return discretized great-circle arcs for each unique undirected edge,
as a vector of `Point3f` arrays. Edges are derived from `cell_nodes`
boundary order, so no `cell_edges` ordering is assumed.
Degenerate edges (coincident endpoints) return a single-point segment.
"""
function edge_segments(g::AbstractManifoldMesh; n_arc_points::Int = 20)
    seen = Set{Tuple{Int, Int}}()
    segments = Vector{Vector{Point3f}}()
    for cid in 1:num_cells(g)
        ns = cell_nodes(g, cid)
        K = length(ns)
        for i in 1:K
            n1, n2 = ns[i], ns[mod1(i + 1, K)]
            key = n1 < n2 ? (n1, n2) : (n2, n1)
            if key ∉ seen
                push!(seen, key)
                p1 = node_coordinates(g, n1)
                p2 = node_coordinates(g, n2)
                push!(segments, slerp(p1, p2, n_arc_points))
            end
        end
    end
    return segments
end

"""
    cell_polygons(g::AbstractManifoldMesh; n_arc_points::Int=20) -> Vector{Vector{Point3f}}

Return closed cell boundaries as discretized great-circle polygons.
Each polygon is a vector of `Point3f` where the first point equals the last.
Edges are derived from consecutive `cell_nodes` in boundary order, so no
`cell_edges` ordering is assumed.
"""
function cell_polygons(g::AbstractManifoldMesh; n_arc_points::Int = 20)
    polygons = Vector{Vector{Point3f}}(undef, num_cells(g))
    for cid in 1:num_cells(g)
        ns = cell_nodes(g, cid)
        K = length(ns)
        polygon = Point3f[]
        for i in 1:K
            n1, n2 = ns[i], ns[mod1(i + 1, K)]
            arc = slerp(node_coordinates(g, n1), node_coordinates(g, n2), n_arc_points)
            if i == 1
                append!(polygon, arc)
            else
                append!(polygon, arc[2:end])
            end
        end
        polygons[cid] = polygon
    end
    return polygons
end
