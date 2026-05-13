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
