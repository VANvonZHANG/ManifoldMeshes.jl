@testset "zero-allocation cached queries" begin
    g = LatLonGrid(lat_edges = collect(-90.0:10.0:90.0),
        lon_edges = collect(0.0:15.0:360.0))

    # Warm up JIT
    cell_volume(g, 1)
    cell_centroid(g, 1)
    node_coordinates(g, 1)

    # cell_volume: pre-cached matrix read — must be zero allocation
    @test (@allocated cell_volume(g, 1)) == 0

    # cell_centroid: pre-cached matrix read — must be zero allocation
    @test (@allocated cell_centroid(g, 1)) == 0

    # node_coordinates: pre-cached matrix read — must be zero allocation
    @test (@allocated node_coordinates(g, 1)) == 0
end

@testset "allocation baseline for on-demand queries" begin
    g = LatLonGrid(lat_edges = collect(-90.0:10.0:90.0),
        lon_edges = collect(0.0:15.0:360.0))

    # Warm up JIT
    edge_length(g, 1)
    edge_midpoint(g, 1)

    # Record baseline allocations (non-zero is acceptable, only for regression detection)
    edge_len_alloc = @allocated edge_length(g, 1)
    edge_mp_alloc = @allocated edge_midpoint(g, 1)

    @testset "edge_length baseline" begin
        # Just records the value; the real check is that it doesn't regress
        @test edge_len_alloc < 10_000  # sanity upper bound (bytes)
    end

    @testset "edge_midpoint baseline" begin
        @test edge_mp_alloc < 10_000
    end
end
