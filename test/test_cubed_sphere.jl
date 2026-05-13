using Manifolds
using ManifoldsBase
using StaticArrays
using LinearAlgebra
using Test

# Include core definitions directly (CubedSphereGrid not yet wired into module)
include(joinpath(@__DIR__, "..", "src", "traits.jl"))
include(joinpath(@__DIR__, "..", "src", "interface.jl"))
include(joinpath(@__DIR__, "..", "src", "sphere", "utils.jl"))
include(joinpath(@__DIR__, "..", "src", "sphere", "cubed_sphere.jl"))

@testset "CubedSphereGrid construction" begin
    g = CubedSphereGrid(n = 2)
    @test g.n == 2
    @test g.R == 1.0
    @test TopologyStyle(g) === IsGrid()
    @test CellTypeStyle(g) === IsUniform{4}()
    @test PatchStyle(g) === MultiPatch{6}()
    @test manifold(g) isa Sphere
    @test num_cells(g) == 6 * 2 * 2
end
