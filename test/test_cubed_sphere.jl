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
    @test num_nodes(g) == 26  # 8 corners + 12 edge-midpoints + 6 face-centers
end

@testset "CubedSphereGrid area conservation" begin
    g = CubedSphereGrid(n = 4)
    total = sum(cell_volume(g, i) for i in 1:num_cells(g))
    @test total ≈ 4π * g.R^2 rtol=1e-10
end

@testset "CubedSphereGrid nodes on sphere" begin
    g = CubedSphereGrid(n = 3, R = 2.0)
    for i in 1:num_nodes(g)
        p = node_coordinates(g, i)
        @test norm(p) ≈ 2.0 atol=1e-12
    end
end

@testset "CubedSphereGrid cell_nodes accessible" begin
    g = CubedSphereGrid(n = 2)
    # Every cell's nodes should be valid
    for cell_id in 1:num_cells(g)
        cn = cell_nodes(g, cell_id)
        for nid in cn
            @test 1 <= nid <= num_nodes(g)
        end
    end
end

@testset "CubedSphereGrid patch queries" begin
    g = CubedSphereGrid(n = 3)
    for face in 1:6
        cells_in_face = count(i -> cell_face(g, i) == face, 1:num_cells(g))
        @test cells_in_face == g.n * g.n
    end
    for i in 1:num_cells(g)
        li, lj = cell_local_2d(g, i)
        @test 1 <= li <= g.n
        @test 1 <= lj <= g.n
    end
end

@testset "CubedSphereGrid equiangular projection" begin
    g_eq = CubedSphereGrid(n = 4, projection = :equiangular)
    total_eq = sum(cell_volume(g_eq, i) for i in 1:num_cells(g_eq))
    @test total_eq ≈ 4π rtol=1e-10
end

@testset "CubedSphereGrid rotation" begin
    rot = SMatrix{3, 3}(0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
    g = CubedSphereGrid(n = 2, rotation = rot)
    total = sum(cell_volume(g, i) for i in 1:num_cells(g))
    @test total ≈ 4π rtol=1e-10
end
