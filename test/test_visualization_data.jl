using ManifoldMeshes
using Manifolds
using StaticArrays
using Test
using GeometryBasics: Point3f
using LinearAlgebra

@testset "Visualization data extraction" begin
    grids = [
        ("LatLonGrid",
            LatLonGrid(lat_edges = collect(-90.0:30.0:90.0), lon_edges = collect(0.0:60.0:360.0))),
        ("CubedSphereGrid", CubedSphereGrid(n = 4)),
        ("ReducedGaussianGrid", ReducedGaussianGrid(nlat = 8)),
        ("HEALPixGrid", HEALPixGrid(nside = 2))
    ]

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

        # Antipodal endpoints
        p3 = SVector(1.0, 0.0, 0.0)
        p4 = SVector(-1.0, 0.0, 0.0)
        pts = ManifoldMeshes.slerp(p3, p4, 3)
        @test length(pts) == 3
        @test all(p -> isapprox(norm(p), 1.0f0, atol = 1e-5), pts)
    end

    @testset "node_points" begin
        g = LatLonGrid(lat_edges = collect(-90.0:30.0:90.0), lon_edges = collect(0.0:60.0:360.0))
        pts = ManifoldMeshes.node_points(g)
        @test pts isa Vector{Point3f}
        @test length(pts) == num_nodes(g)
        @test isapprox(pts[1], Point3f(0, 0, -1), atol = 1e-6)
        @test all(p -> isapprox(norm(p), 1.0f0, atol = 1e-5), pts)
    end

    @testset "edge_segments on $name" for (name, g) in grids
        segs = edge_segments(g; n_arc_points = 10)
        # LatLonGrid has periodic boundary where lon=0 and lon=360 are
        # the same physical line but different node IDs, so edge_segments
        # deduplication sees them as distinct; other grids merge shared nodes
        if name == "LatLonGrid"
            @test length(segs) >= num_edges(g)
        else
            @test length(segs) == num_edges(g)
        end
        for seg in segs
            @test all(!isnan, seg)
            @test length(seg) >= 1
        end
    end

    @testset "cell_polygons on $name" for (name, g) in grids
        polys = cell_polygons(g; n_arc_points = 10)
        @test length(polys) == num_cells(g)
        for poly in polys
            @test all(!isnan, poly)
            @test length(poly) >= 3
            if length(poly) > 1
                @test norm(poly[1] - poly[end]) < 1.0f-3
            end
        end
    end

    @testset "cell_triangles on $name" for (name, g) in grids
        verts, faces = cell_triangles(g)
        n_cells_valid = count(cid -> length(cell_nodes(g, cid)) >= 3, 1:num_cells(g))
        @test length(faces) ==
              sum(length(cell_nodes(g, cid))
        for cid in 1:num_cells(g) if length(cell_nodes(g, cid)) >= 3) * 3
        @test all(f -> 1 <= f <= length(verts), faces)
        @test all(!isnan, verts)

        # Spot-check: centroids should lie on the sphere surface
        R = sqrt(sum(node_coordinates(g, 1) .^ 2))
        offset = 0
        for cid in 1:num_cells(g)
            ns = cell_nodes(g, cid)
            K = length(ns)
            K < 3 && continue
            @test isapprox(norm(verts[offset + 1]), R, atol = 1.0f-4)
            offset += K + 1
        end
    end
end
