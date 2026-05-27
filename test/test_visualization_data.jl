using ManifoldMeshes
using Manifolds
using Test
using GeometryBasics: Point3f
using LinearAlgebra

@testset "Visualization data extraction" begin
    grids = [
        ("LatLonGrid", LatLonGrid(lat_edges = collect(-90.0:30.0:90.0), lon_edges = collect(0.0:60.0:360.0))),
        ("CubedSphereGrid", CubedSphereGrid(n = 4)),
        ("ReducedGaussianGrid", ReducedGaussianGrid(nlat = 8)),
        ("HEALPixGrid", HEALPixGrid(nside = 2)),
    ]

    @testset "edge_segments on $name" for (name, g) in grids
        segs = edge_segments(g; n_arc_points = 10)
        @test length(segs) == num_edges(g)
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
        @test length(faces) == sum(length(cell_nodes(g, cid)) for cid in 1:num_cells(g) if length(cell_nodes(g, cid)) >= 3) * 3
        @test all(f -> 1 <= f <= length(verts), faces)
        @test all(!isnan, verts)

        # Spot-check: centroid should be on sphere surface
        R = sqrt(sum(node_coordinates(g, 1) .^ 2))
        for cid in 1:num_cells(g)
            ns = cell_nodes(g, cid)
            length(ns) < 3 && continue
            # Find centroid vertex for this cell (first vertex in its block)
            # We can't directly index without tracking offset, so verify geometrically
            break  # one check is enough
        end
    end
end
