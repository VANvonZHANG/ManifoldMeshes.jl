using ManifoldMeshes
using LinearAlgebra
using StaticArrays
using Test

@testset "side_of_geodesic" begin
    a = SVector(1.0, 0.0, 0.0)
    b = SVector(0.0, 1.0, 0.0)
    north = SVector(0.0, 0.0, 1.0)
    south = SVector(0.0, 0.0, -1.0)
    mid = SVector(1.0, 1.0, 0.0) / sqrt(2)

    @test side_of_geodesic(north, a, b) == 1
    @test side_of_geodesic(south, a, b) == -1
    @test side_of_geodesic(north, b, a) == -1        # reversing the arc flips the side
    @test side_of_geodesic(a, a, b) == 0             # an endpoint lies on the arc
    @test side_of_geodesic(mid, a, b) == 0           # an interior arc point lies on the arc
    @test side_of_geodesic(2.0 * north, a, b) == 1   # the sign is scale invariant
    @test side_of_geodesic(north, a, a) == 0         # degenerate arc
end

@testset "geodesic_arc_intersection" begin
    a1 = SVector(1.0, 0.0, 0.0)
    b1 = SVector(0.0, 1.0, 0.0)
    a2 = SVector(0.0, 0.0, 1.0)
    b2 = SVector(1.0, 1.0, 0.0) / sqrt(2)
    expected = SVector(1.0, 1.0, 0.0) / sqrt(2)

    x = geodesic_arc_intersection(a1, b1, a2, b2)
    @test x !== nothing
    @test x ≈ expected
    @test geodesic_arc_intersection(a2, b2, a1, b1) ≈ expected   # argument order
    @test geodesic_arc_intersection(b1, a1, b2, a2) ≈ expected   # arc direction

    # the lon = 135° meridian meets the equator at lon = 135°, outside lon ∈ [0°, 90°]
    a3 = SVector(0.0, 0.0, 1.0)
    b3 = SVector(-1.0, 1.0, 0.0) / sqrt(2)
    @test geodesic_arc_intersection(a1, b1, a3, b3) === nothing

    # arcs on one great circle share no isolated crossing point
    a4 = SVector(0.0, 1.0, 0.0)
    b4 = SVector(-1.0, 0.0, 0.0)
    @test geodesic_arc_intersection(a1, b1, a4, b4) === nothing
end

@testset "cell_ring" begin
    g = LatLonGrid(
        lat_edges = [-90.0, -45.0, 45.0, 90.0],
        lon_edges = collect(0.0:90.0:360.0)
    )
    rings = [cell_ring(g, c) for c in 1:num_cells(g)]

    @test all(r -> all(p -> norm(p) ≈ 1.0, r), rings)
    # the middle latitude band is 4 quads; both polar bands collapse to triangles
    @test count(r -> length(r) == 4, rings) == 4
    @test count(r -> length(r) == 3, rings) == 8
end
