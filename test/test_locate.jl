using ManifoldMeshes
using StaticArrays
using Test

@testset "locate helpers" begin
    # _bilinear_weights: corners are unit vectors
    @test ManifoldMeshes._bilinear_weights(0.0, 0.0) == (1.0, 0.0, 0.0, 0.0)   # SW
    @test ManifoldMeshes._bilinear_weights(1.0, 0.0) == (0.0, 1.0, 0.0, 0.0)   # SE
    @test ManifoldMeshes._bilinear_weights(1.0, 1.0) == (0.0, 0.0, 1.0, 0.0)   # NE
    @test ManifoldMeshes._bilinear_weights(0.0, 1.0) == (0.0, 0.0, 0.0, 1.0)   # NW
    # centroid: all 0.25
    @test ManifoldMeshes._bilinear_weights(0.5, 0.5) == (0.25, 0.25, 0.25, 0.25)
    # partition of unity
    for (s, t) in [(0.3, 0.7), (0.9, 0.1), (0.123, 0.456)]
        @test sum(ManifoldMeshes._bilinear_weights(s, t)) ≈ 1.0
    end
end

@testset "coordinate round-trip" begin
    for (lat,
        lon) in [(0.0, 0.0), (45.0, 90.0), (-30.0, 200.0), (89.9, -175.0), (0.0, 359.9)]
        p = ManifoldMeshes._latlon_to_cartesian(lat, lon, 1.0)
        lat2, lon2 = ManifoldMeshes._cartesian_to_latlon(p)
        @test lat2 ≈ lat atol = 1e-9
        @test min(mod(lon2 - lon, 360.0), 360.0 - mod(lon2 - lon, 360.0)) ≈ 0 atol = 1e-9
    end
    # pole: longitude is arbitrary but latitude exact
    p = ManifoldMeshes._latlon_to_cartesian(90.0, 42.0)
    lat2, lon2 = ManifoldMeshes._cartesian_to_latlon(p)
    @test lat2 ≈ 90.0 atol = 1e-9
end

using ManifoldMeshes: LatLonGrid, num_cells, cell_nodes, cell_centroid,
                      locate_cell, interpolation_weights

@testset "LatLon locate_cell" begin
    g = LatLonGrid(lat_edges = [-90.0, -45.0, 0.0, 45.0, 90.0],
        lon_edges = collect(0.0:90.0:360.0))   # nlat=4, nlon=4
    # centroid round-trip
    for cid in [1, 5, 8, 12, 16]
        c = cell_centroid(g, cid)
        lat, lon = ManifoldMeshes._cartesian_to_latlon(c)
        @test locate_cell(g, lat, lon) == cid
    end
    # pole lands in the polar cap band (topmost)
    @test locate_cell(g, 90.0, 0.0) == 13   # (ilat=4, ilon=1): (4-1)*4+1
    @test locate_cell(g, -90.0, 123.0) == 2 # south pole band, ilon=2 (lon=123 in [90,180))
    # longitude wrap
    @test locate_cell(g, 0.0, 0.0) == locate_cell(g, 0.0, 360.0)
    @test locate_cell(g, 0.0, -0.5) == locate_cell(g, 0.0, 359.5)
    # half-open: south edge belongs to the lower cell
    @test locate_cell(g, -45.0, 0.0) == 5   # (ilat=2,ilon=1) -> lat_edges[2]=-45 inclusive
    # bad latitude throws
    @test_throws ArgumentError locate_cell(g, 90.1, 0.0)
end

@testset "LatLon interpolation_weights" begin
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 90.0],
        lon_edges = collect(0.0:90.0:360.0))
    cid = 1   # (ilat=1, ilon=1): lat [-90,0), lon [0,90)
    nodes, w = interpolation_weights(g, cid, -45.0, 45.0)   # centroid -> (0.5,0.5)
    @test nodes == cell_nodes(g, cid)
    @test collect(w) ≈ [0.25, 0.25, 0.25, 0.25]
    # corner SW: (s,t)=(0,0)
    nodes, w = interpolation_weights(g, cid, -90.0, 0.0)
    @test collect(w) ≈ [1.0, 0.0, 0.0, 0.0]
    # reproduces a known bilinear function f(s,t) = 2 + 3s + 5t + 7s*t.
    # cell_nodes ordering is (SW, SE, NE, NW) -> corner (s,t) ((0,0),(1,0),(1,1),(0,1))
    # -> values (2, 5, 17, 7).
    f(s, t) = 2 + 3s + 5t + 7s * t
    corner_vals = (f(0, 0), f(1, 0), f(1, 1), f(0, 1))
    for (s, t) in [(0.25, 0.5), (0.1, 0.9), (0.7, 0.3), (0.5, 0.5)]
        lat = -90.0 + 90.0 * s      # band [-90,0): s in [0,1]
        lon = 90.0 * t              # band [0,90)
        nodes2, w2 = interpolation_weights(g, cid, lat, lon)
        nmax = maximum(nodes2)
        field = zeros(nmax)
        for (nid, v) in zip(nodes2, corner_vals)
            field[nid] = v
        end
        gathered = sum(w2 .* field[collect(nodes2)])
        @test gathered ≈ f(s, t)
    end
end

using ManifoldMeshes: CubedSphereGrid

@testset "CubedSphere locate_cell" begin
    g = CubedSphereGrid(n = 4)
    # centroid round-trip
    for cid in [1, 6 * 16 ÷ 2, 6 * 16]
        c = cell_centroid(g, cid)
        lat, lon = ManifoldMeshes._cartesian_to_latlon(c)
        @test locate_cell(g, lat, lon) == cid
    end
    # rotated grid: locate in the rotated frame still works (centroid round-trip)
    rot = SMatrix{3, 3}(0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
    gr = CubedSphereGrid(n = 3, rotation = rot)
    for cid in [1, 3 * 9 ÷ 2, 6 * 9]
        c = cell_centroid(gr, cid)
        lat, lon = ManifoldMeshes._cartesian_to_latlon(c)
        @test locate_cell(gr, lat, lon) == cid
    end
    # node round-trip: interior nodes lie in a cell that has them as a corner.
    # Face-boundary nodes are ambiguous (non-deduplicated design: same physical
    # point has different IDs per face), so test strictly-interior nodes.
    for nid in [7, 69, 144]
        p = node_coordinates(g, nid)
        lat, lon = ManifoldMeshes._cartesian_to_latlon(p)
        cid = locate_cell(g, lat, lon)
        @test nid ∈ cell_nodes(g, cid)
    end
end

@testset "CubedSphere interpolation_weights" begin
    g = CubedSphereGrid(n = 2)
    cid = 1
    nodes,
    w = interpolation_weights(g, cid,
        ManifoldMeshes._cartesian_to_latlon(cell_centroid(g, cid))...)
    @test nodes == cell_nodes(g, cid)
    @test sum(w) ≈ 1.0
    @test all(>(0), w)
end

using ManifoldMeshes: HEALPixGrid

@testset "HEALPix locate_cell (nested)" begin
    g = HEALPixGrid(nside = 4, ordering = :nested)
    for cid in [1, 12 * 16 ÷ 2, 12 * 16]
        c = cell_centroid(g, cid)
        lat, lon = ManifoldMeshes._cartesian_to_latlon(c)
        @test locate_cell(g, lat, lon) == cid
    end
end

@testset "HEALPix locate_cell (ring)" begin
    g = HEALPixGrid(nside = 4, ordering = :ring)
    for cid in [1, 12 * 16 ÷ 2, 12 * 16]
        c = cell_centroid(g, cid)
        lat, lon = ManifoldMeshes._cartesian_to_latlon(c)
        @test locate_cell(g, lat, lon) == cid
    end
end

@testset "HEALPix interpolation_weights" begin
    g = HEALPixGrid(nside = 2, ordering = :nested)
    cid = 5
    nodes,
    w = interpolation_weights(g, cid,
        ManifoldMeshes._cartesian_to_latlon(cell_centroid(g, cid))...)
    @test nodes == cell_nodes(g, cid)
    @test sum(w) ≈ 1.0
    @test all(>(0), w)
end

using ManifoldMeshes: ReducedGaussianGrid

@testset "ReducedGaussian locate_cell" begin
    g = ReducedGaussianGrid(nlat = 8)
    for cid in [1, num_cells(g) ÷ 2, num_cells(g)]
        c = cell_centroid(g, cid)
        lat, lon = ManifoldMeshes._cartesian_to_latlon(c)
        @test locate_cell(g, lat, lon) == cid
    end
    @test_throws ArgumentError locate_cell(g, 90.1, 0.0)
end

@testset "ReducedGaussian interpolation_weights" begin
    g = ReducedGaussianGrid(nlat = 6)
    cid = num_cells(g) ÷ 3
    nodes,
    w = interpolation_weights(g, cid,
        ManifoldMeshes._cartesian_to_latlon(cell_centroid(g, cid))...)
    @test nodes == cell_nodes(g, cid)
    @test sum(w) ≈ 1.0
    @test all(>(0), w)
end
