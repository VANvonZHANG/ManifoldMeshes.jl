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

@testset "clipping real cells" begin
    # The clipper's unit tests work on hand-built octants and gnomonic squares.
    # Exercise it on the rings the library actually produces, where the
    # degeneracies live (polar quads that collapse to triangles, seam nodes that
    # dedup away, and shared edges that must clip to nothing).
    grids = (
        LatLonGrid(
            lat_edges = collect(-90.0:30.0:90.0),
            lon_edges = collect(0.0:45.0:360.0)
        ),
        CubedSphereGrid(n = 4),
        ReducedGaussianGrid(nlat = 6)
    )
    for g in grids
        label = string(nameof(typeof(g)))
        @testset "self-clip reproduces cell_volume ($label)" begin
            for c in 1:num_cells(g)
                ring = cell_ring(g, c)
                area = spherical_polygon_area(
                    spherical_polygon_intersection(ring, ring), g.R
                )
                @test area ≈ cell_volume(g, c) rtol = 1e-12
            end
        end
        @testset "adjacent cells do not overlap ($label)" begin
            for c in 1:num_cells(g)
                ring_c = cell_ring(g, c)
                for n in cell_cells(g, c)
                    (n > c && n <= num_cells(g)) || continue
                    ring_n = cell_ring(g, n)
                    # exactly zero, in both orders: the shared edge clips away
                    # to a degenerate ring rather than to a sliver of area
                    @test spherical_polygon_area(
                        spherical_polygon_intersection(ring_c, ring_n), g.R
                    ) == 0.0
                    @test spherical_polygon_area(
                        spherical_polygon_intersection(ring_n, ring_c), g.R
                    ) == 0.0
                end
            end
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

    # The other half of the regression: the degenerate guard must not have been
    # bought by breaking the ordinary case. Cross-check a non-degenerate
    # triangle against Girard's theorem — area = angular excess — which is an
    # independent route to the same number, not a second call to the function.
    interior_angle(a, b, c) = begin
        n1 = cross(a, b)
        n2 = cross(a, c)
        atan(norm(cross(n1, n2)), dot(n1, n2))
    end
    girard(a, b, c) = interior_angle(a, b, c) + interior_angle(b, c, a) +
                      interior_angle(c, a, b) - π

    P = normalize(SVector(0.3, -0.4, 0.8))
    Q = normalize(SVector(-0.2, 0.9, 0.1))
    S = normalize(SVector(0.6, 0.2, -0.5))
    @test spherical_triangle_area(1.0, P, Q, S) ≈ girard(P, Q, S) rtol = 1e-12
    @test spherical_triangle_area(3.0, P, Q, S) ≈ 9 * girard(P, Q, S) rtol = 1e-12
end

@testset "spherical_polygon_intersection" begin
    # Tri-rectangular octant A = {x, y, z ≥ 0}; the lon = 45° meridian plane
    # bisects it, so B is exactly half of A and A ∩ B is B.
    A = [SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0)]
    B = [SVector(1.0, 0.0, 0.0), SVector(1.0, 1.0, 0.0) / sqrt(2), SVector(0.0, 0.0, 1.0)]

    @test spherical_polygon_area(B) ≈ π / 4
    @test spherical_polygon_area(spherical_polygon_intersection(A, B)) ≈ π / 4
    @test spherical_polygon_area(spherical_polygon_intersection(B, A)) ≈ π / 4

    # clipping a ring by itself returns the ring
    self = spherical_polygon_intersection(A, A)
    @test spherical_polygon_area(self) ≈ π / 2

    # antipodal octants meet only at the origin: empty intersection.
    # Negating the vertices alone would flip the ring clockwise (winding -1),
    # violating the clipper's counter-clockwise precondition; `reverse` restores
    # it (`cross(-a, -b) == cross(a, b)`, so a clockwise clip ring describes the
    # same half-spaces as A and would return A unchanged).
    anti = reverse([-p for p in A])
    @test isempty(spherical_polygon_intersection(A, anti))
    @test spherical_polygon_area(spherical_polygon_intersection(A, anti)) == 0.0

    # neighbouring octant shares only the x = 0 meridian arc: zero area
    C = [SVector(0.0, 1.0, 0.0), SVector(-1.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0)]
    @test spherical_polygon_area(spherical_polygon_intersection(A, C)) ≈ 0.0 atol = 1e-15

    # result is a valid counter-clockwise ring of unit vectors
    ring = spherical_polygon_intersection(A, B)
    @test all(p -> norm(p) ≈ 1.0, ring)
    @test length(ring) == 3
    @test dot(ring[1], cross(ring[2], ring[3])) > 0   # CCW, as Task 7 consumes it

    # Gnomonic squares: great circles map to straight lines, so a square whose
    # corner pokes past a vertex of the other square is represented exactly.
    gsq(u0,
        u1,
        v0,
        v1) = [normalize(SVector(u, v, 1.0))
               for (u, v) in ((u0, v0), (u1, v0), (u1, v1), (u0, v1))]
    big = gsq(0.0, 1.0, 0.0, 1.0)
    corner = gsq(0.9, 1.1, -0.1, 0.1)

    # The clipper must be symmetric in its arguments. Clipping only against the
    # clip *arc* (rather than its full great circle) dropped the crossing that
    # falls beyond the arc and cut the poked corner off, which made the two
    # orders disagree by 25.6% (3.6102e-3 vs 2.8747e-3).
    @test spherical_polygon_area(spherical_polygon_intersection(big, corner)) ≈
          spherical_polygon_area(spherical_polygon_intersection(corner, big)) rtol = 1e-12

    # Both orders must also reproduce the true overlap area. In the gnomonic
    # chart the overlap is exactly the rectangle u ∈ [0.9, 1.0], v ∈ [0, 0.1],
    # so its area is the independent quadrature
    # ∫∫ du dv / (1 + u² + v²)^(3/2) = 3.804203843254536e-3 (Simpson with the
    # inner integral in u exact, and a 2D midpoint rule, agree to 12 digits).
    overlap = spherical_polygon_area(spherical_polygon_intersection(big, corner))
    @test overlap > 1.7e-3                       # loose sanity bound
    @test overlap ≈ 3.804203843254536e-3 rtol = 1e-9   # the bound that pins the fix
end
