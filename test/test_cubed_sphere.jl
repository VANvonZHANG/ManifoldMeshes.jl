using ManifoldMeshes
using Manifolds
using StaticArrays
using Test

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
