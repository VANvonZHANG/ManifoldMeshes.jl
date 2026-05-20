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

@testset "ReducedGaussianGrid geometry spot checks" begin
    g = ReducedGaussianGrid(nlat = 6)
    # Check 3 representative cells: first, middle, last
    for i in [1, num_cells(g) ÷ 2, num_cells(g)]
        @test abs(norm(cell_centroid(g, i)) - g.R) < 1e-10
        for n in cell_nodes(g, i)
            @test abs(norm(node_coordinates(g, n)) - g.R) < 1e-10
        end
    end
end

@testset "ReducedGaussianGrid edge enumeration" begin
    g = ReducedGaussianGrid(nlat = 4)
    @test num_edges(g) > 0
    # Spot-check 2 cells have valid edge IDs
    for cid in [1, num_cells(g) ÷ 2]
        edges = cell_edges(g, cid)
        @test edges isa NTuple{4, Int}
        @test all(e -> 1 <= e <= num_edges(g), edges)
    end
end

@testset "ReducedGaussianGrid edge sharing" begin
    g = ReducedGaussianGrid(nlat = 4)
    edge_cell_count = zeros(Int, num_edges(g))
    for i in 1:num_cells(g)
        for e in cell_edges(g, i)
            edge_cell_count[e] += 1
        end
    end
    # Non-degenerate edges must be shared by exactly 2 cells
    non_degen = [e for e in 1:num_edges(g) if g._edge_nodes[e][1] != g._edge_nodes[e][2]]
    @test all(e -> edge_cell_count[e] == 2, non_degen)
end

@testset "ReducedGaussianGrid topology spot checks" begin
    g = ReducedGaussianGrid(nlat = 4)
    # Spot-check cell_cells for 3 representative cells
    for i in [1, num_cells(g) ÷ 2, num_cells(g)]
        neighbors = cell_cells(g, i)
        @test neighbors isa NTuple{4, Int}
        @test all(n -> n == 0 || (1 <= n <= num_cells(g)), neighbors)
    end
    # Symmetry check on 3 cells
    for i in [1, num_cells(g) ÷ 2, num_cells(g)]
        for n in cell_cells(g, i)
            n == 0 && continue
            @test i in cell_cells(g, n)
        end
    end
    # Spot-check node_cells for 3 nodes
    for i in [1, num_nodes(g) ÷ 2, num_nodes(g)]
        cells = node_cells(g, i)
        @test cells isa Vector{Int}
        @test all(c -> 1 <= c <= num_cells(g), cells)
    end
end

@testset "ReducedGaussianGrid edge geometry spot checks" begin
    g = ReducedGaussianGrid(nlat = 4)
    # Spot-check 3 edges
    for i in [1, num_edges(g) ÷ 2, num_edges(g)]
        @test edge_length(g, i) >= 0
        mp = edge_midpoint(g, i)
        @test abs(norm(mp) - g.R) < 1e-10
    end
    # Spot-check outward_normal for 3 cells
    for cid in [1, num_cells(g) ÷ 2, num_cells(g)]
        e = cell_edges(g, cid)[1]
        result = edge_outward_normal(g, e, cid)
        @test abs(norm(result.base_point) - g.R) < 1e-10
        @test hasproperty(result, :normal)
    end
end
