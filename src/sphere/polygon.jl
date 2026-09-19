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

# Crossing of the geodesic arc p → q with the *full great circle* through a and
# b. Sutherland–Hodgman clips against the clip edge's line, not the clip arc:
# when the subject pokes past a clip vertex the crossing lies beyond the arc
# and must still be inserted, otherwise the ring cuts the corner off and the
# overlap area is under-counted.
@inline function _edge_crossing(p::SVector{3, Float64}, q::SVector{3, Float64},
        a::SVector{3, Float64}, b::SVector{3, Float64})
    x = cross(cross(p, q), cross(a, b))
    nx = norm(x)
    nx < _GEODESIC_EPS && return nothing
    x = x / nx
    _on_arc(x, p, q) && return x
    x = -x
    return _on_arc(x, p, q) ? x : nothing
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

The vertices are **unit** vectors even when `g.R != 1`, so this ring alone
carries no radius: pass the mesh radius explicitly when measuring it —
`spherical_polygon_area(cell_ring(g, c), g.R)`. Relying on the `R = 1.0`
default silently mis-scales a non-unit sphere's area by a factor `R²`.
"""
function cell_ring(g::AbstractManifoldMesh, cell_id::Int)
    nodes = cell_nodes(g, cell_id)
    ring = SVector{3, Float64}[SVector{3, Float64}(node_coordinates(g, n)) for n in nodes]
    return _dedup_ring(ring)
end

"""
    spherical_polygon_area(ring, R = 1.0) -> Float64

Area of the convex geodesic polygon bounded by `ring` on a sphere of radius
`R`, as a geodesic triangle fan from `ring[1]` evaluated with l'Huilier's
formula.

This is the closed-form evaluation of the Stokes boundary integral
`area(Ω) = ∮_{∂Ω} α` with `dα = ω`: on a constant-curvature surface a geodesic
polygon's boundary integral collapses onto the spherical excess of its
constituent triangles (Girard's theorem). Rings are assumed convex and
counter-clockwise; rings with fewer than three vertices have zero area.

`R` defaults to `1.0` — the unit sphere. Since `cell_ring` always returns unit
vectors, measuring a real mesh cell means passing its radius:
`spherical_polygon_area(cell_ring(g, c), g.R)`.
"""
function spherical_polygon_area(ring::AbstractVector{SVector{3, Float64}},
        R::Real = 1.0)
    n = length(ring)
    n < 3 && return 0.0
    a = normalize(ring[1])
    b = normalize(ring[2])
    unit_area = 0.0
    for k in 3:n
        c = normalize(ring[k])
        unit_area += spherical_triangle_area(1.0, a, b, c)
        b = c
    end
    return Float64(R)^2 * unit_area
end

"""
    spherical_polygon_intersection(subject, clip) -> Vector{SVector{3,Float64}}

Intersection of two convex counter-clockwise geodesic rings, as a ring.

The result is symmetric in its arguments. Zero overlap is reported as an empty
*or degenerate* ring — one with fewer than three vertices, which
`spherical_polygon_area` scores as `0.0`; two rings that merely share a boundary
arc can come back as a two-vertex ring rather than an empty one. Callers must
therefore test the **area**, not `isempty`, to decide whether two cells overlap.

Spherical Sutherland–Hodgman: clip `subject` by each edge of `clip` in turn,
keeping the part of the ring inside the edge's hemisphere
`{q : (a × b) · q ≥ 0}` and inserting the crossing with the edge's great circle
whenever a ring edge leaves or enters (`_edge_crossing` — the full great circle,
not just the arc, so that a subject poking past a clip vertex keeps its corner).
A vertex within `_GEODESIC_EPS` of the clip edge counts as inside, matching the
half-open tie-break of `locate_cell`.

Both rings must be convex and lie within an open hemisphere, which holds for
every mesh cell in this package.
"""
function spherical_polygon_intersection(subject::AbstractVector{SVector{3, Float64}},
        clip::AbstractVector{SVector{3, Float64}})
    length(clip) < 3 && return SVector{3, Float64}[]
    output = collect(subject)
    for i in eachindex(clip)
        a = clip[i]
        b = clip[mod1(i + 1, length(clip))]
        length(output) < 3 && return SVector{3, Float64}[]
        input = output
        output = SVector{3, Float64}[]
        for j in eachindex(input)
            p = input[j]
            q = input[mod1(j + 1, length(input))]
            sp = side_of_geodesic(p, a, b)
            sq = side_of_geodesic(q, a, b)
            if sp >= 0
                if sq >= 0
                    push!(output, q)
                else
                    x = _edge_crossing(p, q, a, b)
                    x === nothing || push!(output, x)
                end
            elseif sq >= 0
                x = _edge_crossing(p, q, a, b)
                x === nothing || push!(output, x)
                push!(output, q)
            end
        end
    end
    return _dedup_ring(output)
end
