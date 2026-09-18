# Shared spherical geometry utilities for S² mesh types

# Angle between two vectors via atan2, numerically stable when they are nearly
# parallel: acos(dot(y, z) / (R * R)) loses half its significant digits there
# (acos(1 - ε) ≈ √(2ε)), which injected spurious area into near-degenerate
# triangles.
@inline _vector_angle(y::SVector{3, Float64}, z::SVector{3, Float64}) = atan(norm(cross(y, z)), dot(y, z))

"""
    spherical_triangle_area(R, A, B, C) -> Float64

Compute the area of a spherical triangle with vertices A, B, C on a sphere
of radius R using l'Huilier's formula.

All vertices must be 3D vectors of length R (i.e., lie on the sphere surface).
"""
function spherical_triangle_area(R::Float64, A::SVector{3, Float64},
        B::SVector{3, Float64}, C::SVector{3, Float64})
    a = _vector_angle(B, C)
    b = _vector_angle(A, C)
    c = _vector_angle(A, B)
    s = (a + b + c) / 2
    tan_half = tan(s / 2) * tan((s - a) / 2) * tan((s - b) / 2) * tan((s - c) / 2)
    E = 4 * atan(sqrt(max(tan_half, 0.0)))
    return R^2 * E
end
