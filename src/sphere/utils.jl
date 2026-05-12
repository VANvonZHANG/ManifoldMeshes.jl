# Shared spherical geometry utilities for S² mesh types

"""
    spherical_triangle_area(R, A, B, C) -> Float64

Compute the area of a spherical triangle with vertices A, B, C on a sphere
of radius R using l'Huilier's formula.

All vertices must be 3D vectors of length R (i.e., lie on the sphere surface).
"""
function spherical_triangle_area(R::Float64, A::SVector{3, Float64},
        B::SVector{3, Float64}, C::SVector{3, Float64})
    a = acos(clamp(dot(B, C) / (R * R), -1, 1))
    b = acos(clamp(dot(A, C) / (R * R), -1, 1))
    c = acos(clamp(dot(A, B) / (R * R), -1, 1))
    s = (a + b + c) / 2
    tan_half = tan(s / 2) * tan((s - a) / 2) * tan((s - b) / 2) * tan((s - c) / 2)
    E = 4 * atan(sqrt(max(tan_half, 0.0)))
    return R^2 * E
end
