# Point location and interpolation-weight primitives.
#
# Per-grid dispatch: each grid implements two private methods
#   _locate_cell(g, lat, lon) -> cell_id
#   _cell_local_coords(g, cell_id, lat, lon) -> (s, t)   # in [0,1)^2
# The shared logic here is coordinate conversion, the bilinear weight formula,
# and the default `interpolation_weights` that composes cell_nodes + weights.

using StaticArrays: SVector
using LinearAlgebra: norm

# ---- Coordinate conversion ----

"""
    _latlon_to_cartesian(lat, lon, R=1.0) -> SVector{3,Float64}

Latitude/longitude in degrees to a Cartesian point on a sphere of radius `R`.
"""
@inline function _latlon_to_cartesian(lat::Real, lon::Real, R::Real = 1.0)
    θ = deg2rad(lat)
    φ = deg2rad(lon)
    sθ, cθ = sincos(θ)
    sφ, cφ = sincos(φ)
    return R * SVector{3, Float64}(cθ * cφ, cθ * sφ, sθ)
end

"""
    _cartesian_to_latlon(p) -> (lat, lon)

Cartesian point (any length) to (latitude, longitude) in degrees.
Longitude normalized to [0, 360).
"""
@inline function _cartesian_to_latlon(p::SVector{3})
    n = norm(p)
    lat = asind(p[3] / n)
    lon = rad2deg(atan(p[2], p[1]))
    if lon < 0
        lon += 360.0
    end
    return (lat, lon)
end

# ---- Bilinear weights ----

"""
    _bilinear_weights(s, t) -> NTuple{4,Float64}

Bilinear weights at local coords `(s, t) ∈ [0,1]²` for a quadrilateral whose
corner node ordering is `(SW, SE, NE, NW)` (matching `cell_nodes`):
`w_SW=(1-s)(1-t)`, `w_SE=s(1-t)`, `w_NE=s*t`, `w_NW=(1-s)*t`.
"""
@inline function _bilinear_weights(s::Real, t::Real)
    return (
        (1 - s) * (1 - t),
        s * (1 - t),
        s * t,
        (1 - s) * t
    )
end

# ---- Public interface ----

"""
    locate_cell(g, lat, lon) -> Int

Return the 1-based cell ID containing the point `(lat, lon)` [degrees].

Latitude must be in `[-90, 90]`; longitude is any real (normalized mod 360).
Tie-break is a half-open `[low, high)` convention per grid type.
"""
function locate_cell(g::AbstractManifoldMesh, lat::Real, lon::Real)
    error("$(typeof(g)) must implement `_locate_cell`")
end

"""
    locate_cell(g, p::SVector{3}) -> Int

3D Cartesian overload: convert `p` to `(lat, lon)` and delegate. `p` need not
be unit length. Intended for mesh-to-mesh workflows where the caller already
holds another mesh's `node_coordinates`.
"""
function locate_cell(g::AbstractManifoldMesh, p::SVector{3})
    lat, lon = _cartesian_to_latlon(p)
    return locate_cell(g, lat, lon)
end

"""
    _cell_local_coords(g, cell_id, lat, lon) -> Tuple{Float64, Float64}

Local coordinates `(s, t) ∈ [0,1]²` of `(lat, lon)` within `cell_id`.
Assumes the point is inside `cell_id` (no check).
"""
function _cell_local_coords(g::AbstractManifoldMesh, cell_id::Int, lat::Real, lon::Real)
    error("$(typeof(g)) must implement `_cell_local_coords`")
end

"""
    interpolation_weights(g, cell_id, lat, lon) -> (nodes, weights)

Bilinear interpolation weights for `(lat, lon)` within `cell_id`.

Returns `(nodes, weights)` where `nodes` matches `cell_nodes(g, cell_id)` in
order and `weights` is the NodeLoc bilinear weight per node. Both are
`NTuple{4}` (zero allocation). The point is assumed to be inside `cell_id`.
"""
function interpolation_weights(g::AbstractManifoldMesh, cell_id::Int, lat::Real, lon::Real)
    nodes = cell_nodes(g, cell_id)
    s, t = _cell_local_coords(g, cell_id, lat, lon)
    return (nodes, _bilinear_weights(s, t))
end
