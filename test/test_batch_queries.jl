using ManifoldMeshes
using Test
using StaticArrays
using LinearAlgebra

@testset "Batch queries" begin
    # HEALPixGrid
    g = HEALPixGrid(nside=4, R=1.0)
    vols = all_cell_volumes(g)
    @test length(vols) == num_cells(g)
    @test sum(vols) ≈ 4π * g.R^2 rtol = 1e-2
    @test vols === g.cell_volumes  # zero-copy

    nodes = all_node_coordinates(g)
    @test length(nodes) == num_nodes(g)
    @test eltype(nodes) == SVector{3, Float64}
    @test nodes === g.nodes  # zero-copy

    # LatLonGrid
    g2 = LatLonGrid(lat_edges=[-90.0, -30.0, 30.0, 90.0], lon_edges=collect(0.0:30.0:360.0), R=1.0)
    vols2 = all_cell_volumes(g2)
    @test length(vols2) == num_cells(g2)
    @test sum(vols2) ≈ 4π * g2.R^2 rtol = 1e-10

    nodes2 = all_node_coordinates(g2)
    @test length(nodes2) == num_nodes(g2)

    # CubedSphereGrid
    g3 = CubedSphereGrid(n=4, R=1.0)
    vols3 = all_cell_volumes(g3)
    @test length(vols3) == num_cells(g3)
    @test sum(vols3) ≈ 4π * g3.R^2 rtol = 1e-10

    # ReducedGaussianGrid
    g4 = ReducedGaussianGrid(nlat=8, R=1.0)
    vols4 = all_cell_volumes(g4)
    @test length(vols4) == num_cells(g4)
    @test sum(vols4) ≈ 4π * g4.R^2 rtol = 1e-8

    # Test all_cell_centroids
    cents = all_cell_centroids(g)
    @test length(cents) == num_cells(g)
    @test eltype(cents) == SVector{3, Float64}
    @test all(c -> norm(c) ≈ g.R, cents)

    # Test all_edge_lengths
    edge_lens = all_edge_lengths(g)
    @test length(edge_lens) == num_edges(g)
    @test all(l -> l >= 0, edge_lens)

    # Test on other grid types
    g2 = LatLonGrid(lat_edges=[-90.0, -30.0, 30.0, 90.0], lon_edges=collect(0.0:30.0:360.0), R=1.0)
    @test length(all_cell_centroids(g2)) == num_cells(g2)
    @test length(all_edge_lengths(g2)) == num_edges(g2)

    g3 = CubedSphereGrid(n=4, R=1.0)
    @test length(all_cell_centroids(g3)) == num_cells(g3)
    @test length(all_edge_lengths(g3)) == num_edges(g3)

    g4 = ReducedGaussianGrid(nlat=8, R=1.0)
    @test length(all_cell_centroids(g4)) == num_cells(g4)
    @test length(all_edge_lengths(g4)) == num_edges(g4)
end
