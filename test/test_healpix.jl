using Manifolds
using ManifoldsBase
using StaticArrays
using LinearAlgebra
using Test

# Include core definitions directly (HEALPixGrid not yet wired into module)
include(joinpath(@__DIR__, "..", "src", "traits.jl"))
include(joinpath(@__DIR__, "..", "src", "interface.jl"))
include(joinpath(@__DIR__, "..", "src", "sphere", "utils.jl"))
include(joinpath(@__DIR__, "..", "src", "sphere", "healpix.jl"))

@testset "HEALPixGrid construction" begin
    g = HEALPixGrid(nside = 2)
    @test g.nside == 2
    @test TopologyStyle(g) === IsSemiGrid()
    @test CellTypeStyle(g) === IsUniform{4}()
    @test PatchStyle(g) === NoPatch()
    @test manifold(g) isa Sphere
    @test num_cells(g) == 12 * 2 * 2
end

@testset "HEALPixGrid cell count invariant" begin
    for nside in [1, 2, 4, 8]
        g = HEALPixGrid(nside = nside)
        @test num_cells(g) == 12 * nside^2
    end
end

@testset "HEALPixGrid area conservation" begin
    g = HEALPixGrid(nside = 4)
    total = sum(cell_volume(g, i) for i in 1:num_cells(g))
    @test total ≈ 4π * g.R^2 rtol = 1e-10
end
