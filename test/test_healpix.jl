using ManifoldMeshes
using Manifolds
using ManifoldsBase
using StaticArrays
using LinearAlgebra
using Test

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

@testset "HEALPixGrid centroids on sphere" begin
    g = HEALPixGrid(nside = 4)
    for i in 1:num_cells(g)
        c = cell_centroid(g, i)
        @test abs(norm(c) - g.R) < 1e-10
    end
end

@testset "HEALPixGrid cell_nodes" begin
    g = HEALPixGrid(nside = 2)
    for i in 1:num_cells(g)
        nodes = cell_nodes(g, i)
        @test nodes isa NTuple{4, Int}
        for n in nodes
            p = node_coordinates(g, n)
            @test abs(norm(p) - g.R) < 1e-10
        end
    end
end

@testset "HEALPixGrid nodes are unique vertices" begin
    g = HEALPixGrid(nside = 2)
    coords = [node_coordinates(g, i) for i in 1:num_nodes(g)]
    rounded = [round.(c, digits = 10) for c in coords]
    @test length(unique(rounded)) == num_nodes(g)
end

@testset "HEALPixGrid num_edges" begin
    g = HEALPixGrid(nside = 2)
    @test num_edges(g) > 0
end

@testset "HEALPixGrid cell_edges" begin
    g = HEALPixGrid(nside = 2)
    for i in 1:num_cells(g)
        edges = cell_edges(g, i)
        @test edges isa NTuple{4, Int}
    end
end

@testset "HEALPixGrid cell_cells" begin
    g = HEALPixGrid(nside = 2)
    for i in 1:num_cells(g)
        neighbors = cell_cells(g, i)
        @test neighbors isa NTuple{4, Int}
    end
end

@testset "HEALPixGrid node_cells" begin
    g = HEALPixGrid(nside = 2)
    for i in 1:num_nodes(g)
        cells = node_cells(g, i)
        @test cells isa Vector{Int}
    end
end

@testset "HEALPixGrid edge geometry" begin
    g = HEALPixGrid(nside = 2)
    for i in 1:num_edges(g)
        @test edge_length(g, i) >= 0
        mp = edge_midpoint(g, i)
        @test abs(norm(mp) - g.R) < 1e-10
    end

    for i in 1:num_cells(g)
        for e in cell_edges(g, i)
            result = edge_outward_normal(g, e, i)
            @test abs(norm(result.base_point) - g.R) < 1e-10
        end
    end
end
