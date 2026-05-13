using ManifoldMeshes
using Manifolds
using ManifoldsBase
using StaticArrays
using LinearAlgebra
using Test

@testset "ReducedGaussianGrid construction" begin
    g = ReducedGaussianGrid(nlat = 4)
    @test g.nlat == 4
    @test TopologyStyle(g) === IsSemiGrid()
    @test CellTypeStyle(g) === IsUniform{4}()
    @test PatchStyle(g) === NoPatch()
    @test manifold(g) isa Sphere
    @test num_cells(g) > 0
end

@testset "ReducedGaussianGrid area conservation" begin
    g = ReducedGaussianGrid(nlat = 8)
    total = sum(cell_volume(g, i) for i in 1:num_cells(g))
    @test total ≈ 4π * g.R^2 rtol = 1e-8
end

@testset "ReducedGaussianGrid centroids on sphere" begin
    g = ReducedGaussianGrid(nlat = 6)
    for i in 1:num_cells(g)
        c = cell_centroid(g, i)
        @test abs(norm(c) - g.R) < 1e-10
    end
end

@testset "ReducedGaussianGrid cell_nodes" begin
    g = ReducedGaussianGrid(nlat = 4)

    for i in 1:num_cells(g)
        nodes = cell_nodes(g, i)
        @test nodes isa NTuple{4, Int}
        for node_id in nodes
            p = node_coordinates(g, node_id)
            @test abs(norm(p) - g.R) < 1e-10
        end
    end
end

@testset "ReducedGaussianGrid num_edges" begin
    g = ReducedGaussianGrid(nlat = 4)
    @test num_edges(g) > 0
end

@testset "ReducedGaussianGrid cell_edges" begin
    g = ReducedGaussianGrid(nlat = 4)
    for i in 1:num_cells(g)
        edges = cell_edges(g, i)
        @test edges isa NTuple{4, Int}
        @test all(e -> 1 <= e <= num_edges(g), edges)
    end
end

@testset "ReducedGaussianGrid edge consistency" begin
    g = ReducedGaussianGrid(nlat = 4)
    edge_cell_count = zeros(Int, num_edges(g))
    for i in 1:num_cells(g)
        for e in cell_edges(g, i)
            edge_cell_count[e] += 1
        end
    end
    # Non-degenerate edges (distinct endpoints) must be shared by exactly 2 cells.
    # Self-loop edges (collapsed endpoints) occur at poles and at latitude-band
    # transitions where node counts differ; they are valid but have different sharing.
    for e in 1:num_edges(g)
        n1, n2 = g._edge_nodes[e]
        if n1 != n2
            @test edge_cell_count[e] == 2
        end
    end
end
