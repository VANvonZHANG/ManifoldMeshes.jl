using ManifoldMeshes
using Test
using StaticArrays

@testset "node_cells correctness" begin
    # --- LatLonGrid: interior node should have 4 adjacent cells ---
    g = LatLonGrid(
        lat_edges = [-90.0, -30.0, 30.0, 90.0],
        lon_edges = collect(0.0:30.0:360.0),
        R = 1.0,
    )
    # Node at intersection of lat=-30, lon=30 (interior)
    # ilat=2, ilon=2 -> node_id = (2-1)*(12+1) + 2 = 15
    node_id = 15
    cells = node_cells(g, node_id)
    @test length(cells) == 4
    # Cells: SW=(1,1), SE=(1,2), NW=(2,1), NE=(2,2)
    @test sort(collect(cells)) == [1, 2, 13, 14]

    # --- CubedSphereGrid: every node should belong to at least 1 cell ---
    g2 = CubedSphereGrid(n = 4, R = 1.0)
    for node_id in 1:num_nodes(g2)
        cells = node_cells(g2, node_id)
        @test length(cells) >= 1
        @test all(c -> 1 <= c <= num_cells(g2), cells)
    end

    # --- HEALPixGrid: every node should belong to at least 1 cell ---
    g3 = HEALPixGrid(nside = 4, R = 1.0)
    for node_id in 1:num_nodes(g3)
        cells = node_cells(g3, node_id)
        @test length(cells) >= 1
        @test all(c -> 1 <= c <= num_cells(g3), cells)
    end

    # --- ReducedGaussianGrid: every node should belong to at least 1 cell ---
    g4 = ReducedGaussianGrid(nlat = 8, R = 1.0)
    for node_id in 1:num_nodes(g4)
        cells = node_cells(g4, node_id)
        @test length(cells) >= 1
        @test all(c -> 1 <= c <= num_cells(g4), cells)
    end

    # Consistency: node_cells result must be stable across repeated calls
    g5 = HEALPixGrid(nside = 2, R = 1.0)
    for node_id in 1:min(20, num_nodes(g5))
        c1 = collect(node_cells(g5, node_id))
        c2 = collect(node_cells(g5, node_id))
        @test c1 == c2
    end
end
