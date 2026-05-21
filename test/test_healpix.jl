using ManifoldMeshes
using Manifolds
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
    @test total ≈ 4π * g.R^2 rtol = 1e-2
end

@testset "HEALPixGrid geometry spot checks" begin
    g = HEALPixGrid(nside = 4)
    # Spot-check 3 cells: first, middle, last
    for i in [1, num_cells(g) ÷ 2, num_cells(g)]
        @test abs(norm(cell_centroid(g, i)) - g.R) < 1e-10
        nodes = cell_nodes(g, i)
        @test nodes isa NTuple{4, Int}
        for n in nodes
            @test abs(norm(node_coordinates(g, n)) - g.R) < 1e-10
        end
    end
end

@testset "HEALPixGrid nodes are unique vertices" begin
    g = HEALPixGrid(nside = 2)
    coords = [node_coordinates(g, i) for i in 1:num_nodes(g)]
    rounded = [round.(c, digits = 10) for c in coords]
    @test length(unique(rounded)) == num_nodes(g)
end

@testset "HEALPixGrid topology spot checks" begin
    g = HEALPixGrid(nside = 2)
    @test num_edges(g) > 0
    # Spot-check 2 cells
    for cid in [1, num_cells(g) ÷ 2]
        @test cell_edges(g, cid) isa NTuple{4, Int}
        @test cell_cells(g, cid) isa NTuple{4, Int}
    end
    # Spot-check 3 nodes
    for i in [1, num_nodes(g) ÷ 2, num_nodes(g)]
        @test node_cells(g, i) isa Vector{Int}
    end
end

@testset "HEALPixGrid edge geometry spot checks" begin
    g = HEALPixGrid(nside = 2)
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
    end
end

@testset "HEALPixGrid nested ordering" begin
    for nside in [1, 2, 4]
        g_ring = HEALPixGrid(nside = nside, ordering = :ring)
        g_nested = HEALPixGrid(nside = nside, ordering = :nested)

        @test g_ring.nside == g_nested.nside
        @test num_cells(g_ring) == num_cells(g_nested)
        @test num_nodes(g_ring) == num_nodes(g_nested)

        total_ring = sum(cell_volume(g_ring, i) for i in 1:num_cells(g_ring))
        total_nested = sum(cell_volume(g_nested, i) for i in 1:num_cells(g_nested))
        @test total_ring ≈ total_nested rtol = 1e-12

        cents_ring = Set(round.(cell_centroid(g_ring, i), digits = 10)
        for i in 1:num_cells(g_ring))
        cents_nested = Set(round.(cell_centroid(g_nested, i), digits = 10)
        for i in 1:num_cells(g_nested))
        @test cents_ring == cents_nested

        vols_ring = sort([cell_volume(g_ring, i) for i in 1:num_cells(g_ring)])
        vols_nested = sort([cell_volume(g_nested, i) for i in 1:num_cells(g_nested)])
        @test vols_ring ≈ vols_nested rtol = 1e-12
    end
end

@testset "HEALPixGrid nested connectivity consistency" begin
    g_ring = HEALPixGrid(nside = 2, ordering = :ring)
    g_nested = HEALPixGrid(nside = 2, ordering = :nested)

    ring_zeros = count(n == 0 for i in 1:num_cells(g_ring) for n in cell_cells(g_ring, i))
    nested_zeros = count(n == 0 for i in 1:num_cells(g_nested)
    for n in cell_cells(g_nested, i))
    @test ring_zeros == nested_zeros

    function count_asymmetric(g)
        count = 0
        for i in 1:num_cells(g)
            for n in cell_cells(g, i)
                if n > 0 && !(i in cell_cells(g, n))
                    count += 1
                end
            end
        end
        return count
    end
    @test count_asymmetric(g_ring) == count_asymmetric(g_nested)

    ring_edge_count = zeros(Int, num_edges(g_ring))
    for c in 1:num_cells(g_ring)
        for e in cell_edges(g_ring, c)
            ring_edge_count[e] += 1
        end
    end
    nested_edge_count = zeros(Int, num_edges(g_nested))
    for c in 1:num_cells(g_nested)
        for e in cell_edges(g_nested, c)
            nested_edge_count[e] += 1
        end
    end
    @test sort(ring_edge_count) == sort(nested_edge_count)
end
