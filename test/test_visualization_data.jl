using ManifoldMeshes
using Manifolds
using StaticArrays
using Test
using GeometryBasics

@testset "Visualization data extraction" begin
    grid = LatLonGrid(lat_edges = Float64.(collect(-90:30:90)), lon_edges = Float64.(collect(0:60:360)))

    @testset "slerp" begin
        # Coincident points -> single point
        p = SVector(0.0, 0.0, 1.0)
        pts = ManifoldMeshes.slerp(p, p, 5)
        @test length(pts) == 1
        @test isapprox(pts[1], Point3f(0, 0, 1))

        # Orthogonal points (equator to north pole)
        p1 = SVector(1.0, 0.0, 0.0)
        p2 = SVector(0.0, 0.0, 1.0)
        pts = ManifoldMeshes.slerp(p1, p2, 3)
        @test length(pts) == 3
        @test isapprox(pts[1], Point3f(Float32.(p1)), atol = 1e-6)
        @test isapprox(pts[3], Point3f(Float32.(p2)), atol = 1e-6)
        # All points on unit sphere
        for pt in pts
            @test isapprox(norm(pt), 1.0f0, atol = 1e-6)
        end

        # All points lie on the same great circle plane (normal = cross(p1, p2))
        normal = cross(p1, p2)
        for pt in pts
            v = SVector(Float64(pt[1]), Float64(pt[2]), Float64(pt[3]))
            @test isapprox(abs(dot(v, normal)), 0.0, atol = 1e-6)
        end
    end

    @testset "node_points" begin
        pts = ManifoldMeshes.node_points(grid)
        @test pts isa Vector{Point3f}
        @test length(pts) == num_nodes(grid)  # 7 × 7 = 49
        # First node at south pole (0, 0, -R) with R=1
        @test isapprox(pts[1][1], 0.0, atol = 1e-10)
        @test isapprox(pts[1][2], 0.0, atol = 1e-10)
        @test isapprox(pts[1][3], -1.0, atol = 1e-10)
        # All points on unit sphere
        for p in pts
            @test isapprox(norm(p), 1.0f0, atol = 1e-6)
        end
    end

    @testset "edge_segments" begin
        segs = ManifoldMeshes.edge_segments(grid, n_arc_points = 5)
        @test segs isa Vector{Vector{Point3f}}
        @test length(segs) == num_edges(grid)
        for seg in segs
            @test length(seg) == 5
        end
        # Non-degenerate edges: all points on unit sphere
        for seg in segs
            if length(unique(seg)) > 1
                for p in seg
                    @test isapprox(norm(p), 1.0f0, atol = 1e-6)
                end
            end
        end
    end

    @testset "cell_polygons" begin
        polys = ManifoldMeshes.cell_polygons(grid, n_arc_points = 5)
        @test polys isa Vector{Vector{Point3f}}
        @test length(polys) == num_cells(grid)
        # Each polygon is closed (first == last)
        for poly in polys
            @test isapprox(poly[1], poly[end], atol = 1e-6)
        end
        # Non-degenerate cells: all points on unit sphere
        for poly in polys
            if length(poly) > 2
                for p in poly
                    @test isapprox(norm(p), 1.0f0, atol = 1e-6)
                end
            end
        end
    end

    @testset "polar degeneracy — no NaN" begin
        segs = ManifoldMeshes.edge_segments(grid, n_arc_points = 10)
        for seg in segs
            for p in seg
                @test !isnan(p[1]) && !isnan(p[2]) && !isnan(p[3])
            end
        end
    end

    @testset "periodic boundary — seam edges" begin
        segs = ManifoldMeshes.edge_segments(grid, n_arc_points = 5)
        for seg in segs
            if length(unique(seg)) > 1
                for p in seg
                    @test isapprox(norm(p), 1.0f0, atol = 1e-6)
                end
            end
        end
    end

    @testset "cell_polygons — no NaN" begin
        polys = ManifoldMeshes.cell_polygons(grid, n_arc_points = 10)
        for poly in polys
            for p in poly
                @test !isnan(p[1]) && !isnan(p[2]) && !isnan(p[3])
            end
        end
    end
end
