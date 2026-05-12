using ManifoldMeshes
using Manifolds
using StaticArrays
using Test
using GeometryBasics
using LinearAlgebra

@testset "Visualization data extraction" begin
    grid = LatLonGrid(lat_edges = Float64.(collect(-90:30:90)), lon_edges = Float64.(collect(0:60:360)))

    @testset "slerp" begin
        p = SVector(0.0, 0.0, 1.0)
        pts = ManifoldMeshes.slerp(p, p, 5)
        @test length(pts) == 1

        p1 = SVector(1.0, 0.0, 0.0)
        p2 = SVector(0.0, 0.0, 1.0)
        pts = ManifoldMeshes.slerp(p1, p2, 3)
        @test length(pts) == 3
        @test isapprox(pts[1], Point3f(Float32.(p1)), atol = 1e-6)
        @test isapprox(pts[3], Point3f(Float32.(p2)), atol = 1e-6)
        @test isapprox(norm(pts[2]), 1.0f0, atol = 1e-6)
    end

    @testset "node_points" begin
        pts = ManifoldMeshes.node_points(grid)
        @test pts isa Vector{Point3f}
        @test length(pts) == num_nodes(grid)
        # Spot check south pole and one equatorial node
        @test isapprox(pts[1], Point3f(0, 0, -1), atol = 1e-6)
        eq_idx = ManifoldMeshes._node_linear_index(grid, 3, 1)
        @test isapprox(norm(pts[eq_idx]), 1.0f0, atol = 1e-6)
    end

    @testset "edge_segments" begin
        segs = ManifoldMeshes.edge_segments(grid, n_arc_points = 5)
        @test segs isa Vector{Vector{Point3f}}
        @test length(segs) == num_edges(grid)
        # Spot check one non-degenerate edge
        @test isapprox(norm(segs[7][3]), 1.0f0, atol = 1e-6)
    end

    @testset "cell_polygons" begin
        polys = ManifoldMeshes.cell_polygons(grid, n_arc_points = 5)
        @test polys isa Vector{Vector{Point3f}}
        @test length(polys) == num_cells(grid)
        # Spot check one polygon is closed
        @test isapprox(polys[3][1], polys[3][end], atol = 1e-6)
    end

    @testset "no NaN in output" begin
        segs = ManifoldMeshes.edge_segments(grid, n_arc_points = 10)
        polys = ManifoldMeshes.cell_polygons(grid, n_arc_points = 10)
        # Spot check a few samples instead of all
        @test all(!isnan, segs[1][1])
        @test all(!isnan, polys[1][1])
    end
end
