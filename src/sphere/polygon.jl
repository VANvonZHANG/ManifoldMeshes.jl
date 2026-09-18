# Spherical polygon geometry: geodesic predicates and ring utilities shared by
# every mesh type.
#
# Ring convention: `Vector{SVector{3,Float64}}` of unit vectors, consecutive
# vertices, counter-clockwise seen from outside the sphere, so
# `side_of_geodesic(q, a, b) > 0` means `q` lies inside the ring containing the
# arc `a → b`. Cells are additionally convex.

using StaticArrays: SVector
using LinearAlgebra: cross, dot, norm, normalize

# Angular tolerance of `side_of_geodesic` (which normalizes its normal, so this
# is radians) and 3D distance tolerance elsewhere in this file.
const _GEODESIC_EPS = 1e-12

"""
    side_of_geodesic(q, a, b) -> Int

Which side of the oriented geodesic arc `a → b` the direction `q` lies on:
`+1` inside the counter-clockwise ring that contains the arc, `-1` outside,
`0` on the arc itself (within `_GEODESIC_EPS`).

The predicate is `sign((a × b) · q)` — the signed volume of the tetrahedron
`(0, a, b, q)`, a coordinate-free orientation test. Intrinsically it is the
"which way does `log_a q` point relative to `a → b`" test, which is what a
general manifold implementation would evaluate in the logarithmic chart.
"""
function side_of_geodesic(q::SVector{3, Float64}, a::SVector{3, Float64},
        b::SVector{3, Float64})
    n = cross(a, b)
    nn = norm(n)
    nn < _GEODESIC_EPS && return 0
    s = dot(n, q) / nn
    abs(s) < _GEODESIC_EPS && return 0
    return s > 0 ? 1 : -1
end

"""
    geodesic_arc_intersection(a1, b1, a2, b2) -> Union{SVector{3,Float64},Nothing}

Crossing point of the geodesic arcs `a1 → b1` and `a2 → b2` (each shorter than
a semicircle), or `nothing` when they do not cross.

The great circles meet at `±normalize((a1 × b1) × (a2 × b2))`; the root that
lies on both arcs is returned. Collinear great circles report `nothing`: a
shared boundary contributes no area.
"""
function geodesic_arc_intersection(a1::SVector{3, Float64}, b1::SVector{3, Float64},
        a2::SVector{3, Float64}, b2::SVector{3, Float64})
    x = cross(cross(a1, b1), cross(a2, b2))
    nx = norm(x)
    nx < _GEODESIC_EPS && return nothing
    x = x / nx
    _on_arc(x, a1, b1) && _on_arc(x, a2, b2) && return x
    x = -x
    _on_arc(x, a1, b1) && _on_arc(x, a2, b2) && return x
    return nothing
end

# `p` lies on the shorter arc `a → b` iff it is neither before `a` nor after
# `b` along the arc orientation, i.e. both `cross(a, p)` and `cross(p, b)`
# point along `cross(a, b)`.
@inline function _on_arc(p::SVector{3, Float64}, a::SVector{3, Float64},
        b::SVector{3, Float64})
    n = cross(a, b)
    tol = _GEODESIC_EPS * norm(n)
    return dot(cross(a, p), n) >= -tol && dot(cross(p, b), n) >= -tol
end

"""
    _dedup_ring(ring) -> Vector{SVector{3,Float64}}

Normalize the vertices of `ring` and drop consecutive duplicates (including
the wrap-around pair).
"""
function _dedup_ring(ring::AbstractVector{SVector{3, Float64}})
    out = SVector{3, Float64}[]
    for p in ring
        q = normalize(p)
        (isempty(out) || norm(out[end] - q) > _GEODESIC_EPS) && push!(out, q)
    end
    if length(out) > 1 && norm(out[1] - out[end]) <= _GEODESIC_EPS
        pop!(out)
    end
    return out
end

"""
    cell_ring(g::AbstractManifoldMesh, cell_id::Int) -> Vector{SVector{3,Float64}}

Boundary of cell `cell_id` as a unit-vector ring in `cell_nodes` order, with
consecutive duplicate vertices removed: `LatLonGrid` polar cells repeat the
pole node and seam cells repeat the `lon = 0°`/`lon = 360°` node, which would
otherwise enter clipping as zero-length edges.

Cells are convex and counter-clockwise seen from outside the sphere.
"""
function cell_ring(g::AbstractManifoldMesh, cell_id::Int)
    nodes = cell_nodes(g, cell_id)
    ring = SVector{3, Float64}[SVector{3, Float64}(node_coordinates(g, n)) for n in nodes]
    return _dedup_ring(ring)
end
