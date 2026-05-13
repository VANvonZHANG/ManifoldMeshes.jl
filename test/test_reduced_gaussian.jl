using Manifolds
using ManifoldsBase
using StaticArrays
using LinearAlgebra
using Test

# Include core definitions directly (ReducedGaussianGrid not yet wired into module)
include(joinpath(@__DIR__, "..", "src", "traits.jl"))
include(joinpath(@__DIR__, "..", "src", "interface.jl"))
include(joinpath(@__DIR__, "..", "src", "sphere", "utils.jl"))
include(joinpath(@__DIR__, "..", "src", "sphere", "reduced_gaussian.jl"))

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
