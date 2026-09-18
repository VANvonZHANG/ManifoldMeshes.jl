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

@testset "spherical_polygon_area matches cell_volume" begin
    face_nodes = [1 2 5; 2 3 5; 3 4 5; 4 1 5; 2 1 6; 3 2 6; 4 3 6; 1 4 6]
    grids = [
        LatLonGrid(
            lat_edges = collect(-90.0:30.0:90.0),
            lon_edges = collect(0.0:45.0:360.0)
        ),
        CubedSphereGrid(n = 4),
        ReducedGaussianGrid(nlat = 6),
        HEALPixGrid(nside = 2),
        UnstructuredMesh(
            [0.0, 90.0, 180.0, 270.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 90.0, -90.0],
            face_nodes;
            start_index = 1
        )
    ]

    for g in grids
        for c in 1:num_cells(g)
            area = spherical_polygon_area(cell_ring(g, c), g.R)
            @test area ≈ cell_volume(g, c) rtol = 1e-10
        end
    end
end

@testset "spherical_polygon_area degenerate rings" begin
    @test spherical_polygon_area(SVector{3, Float64}[]) == 0.0
    p = SVector(1.0, 0.0, 0.0)
    @test spherical_polygon_area([p, p]) == 0.0
    # a tri-rectangular octant triangle has spherical excess π/2
    tri = [SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0)]
    @test spherical_polygon_area(tri) ≈ π / 2
    @test spherical_polygon_area(tri, 2.0) ≈ 4 * π / 2
end

@testset "spherical_triangle_area degenerate triangles" begin
    # nodes 2 and 7 of ReducedGaussianGrid(nlat = 6), the quad of cell 5
    A = SVector(0.74535599249992979, 0.0, -0.66666666666666674)
    B = SVector(0.66666666666666674, 0.66666666666666663, -0.33333333333333337)
    @test norm(A) != 1.0     # one ulp short of the sphere: acos(dot(A, A)) = 1.5e-8
    @test spherical_triangle_area(1.0, A, A, B) == 0.0
    # one ulp apart: still degenerate to well below any plausible tolerance,
    # not the ~6e-9 the acos formulation produced
    @test spherical_triangle_area(1.0, A, nextfloat.(A), B) < 1e-15
end
