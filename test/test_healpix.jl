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
    @test total ≈ 4π * g.R^2 rtol = 1e-3
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

@testset "HEALPixGrid nested ordering" begin
    for nside in [1, 2, 4]
        g_ring = HEALPixGrid(nside = nside, ordering = :ring)
        g_nested = HEALPixGrid(nside = nside, ordering = :nested)

        @test g_ring.nside == g_nested.nside
        @test num_cells(g_ring) == num_cells(g_nested)
        @test num_nodes(g_ring) == num_nodes(g_nested)

        # Total area should be identical
        total_ring = sum(cell_volume(g_ring, i) for i in 1:num_cells(g_ring))
        total_nested = sum(cell_volume(g_nested, i) for i in 1:num_cells(g_nested))
        @test total_ring ≈ total_nested rtol = 1e-12

        # Centroids should be the same set (just reordered)
        cents_ring = Set(round.(cell_centroid(g_ring, i), digits = 10) for i in 1:num_cells(g_ring))
        cents_nested = Set(round.(cell_centroid(g_nested, i), digits = 10) for i in 1:num_cells(g_nested))
        @test cents_ring == cents_nested

        # Cell volumes should be the same multiset
        vols_ring = sort([cell_volume(g_ring, i) for i in 1:num_cells(g_ring)])
        vols_nested = sort([cell_volume(g_nested, i) for i in 1:num_cells(g_nested)])
        @test vols_ring ≈ vols_nested rtol = 1e-12
    end
end

@testset "HEALPixGrid nested ordering connectivity consistency" begin
    g_ring = HEALPixGrid(nside = 2, ordering = :ring)
    g_nested = HEALPixGrid(nside = 2, ordering = :nested)

    # The nested grid should have the same connectivity pattern as the ring grid,
    # just with renumbered cells.

    # Count zero neighbors in both grids
    ring_zeros = count(n == 0 for i in 1:num_cells(g_ring) for n in cell_cells(g_ring, i))
    nested_zeros = count(n == 0 for i in 1:num_cells(g_nested) for n in cell_cells(g_nested, i))
    @test ring_zeros == nested_zeros

    # Count asymmetric neighbor relations in both grids
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

    # Edge sharing pattern should be the same
    ring_edge_count = zeros(Int, num_edges(g_ring))
    for i in 1:num_cells(g_ring)
        for e in cell_edges(g_ring, i)
            ring_edge_count[e] += 1
        end
    end
    nested_edge_count = zeros(Int, num_edges(g_nested))
    for i in 1:num_cells(g_nested)
        for e in cell_edges(g_nested, i)
            nested_edge_count[e] += 1
        end
    end
    @test sort(ring_edge_count) == sort(nested_edge_count)
end

@testset "HEALPixGrid nested ordering cell_nodes valid" begin
    g = HEALPixGrid(nside = 2, ordering = :nested)
    for i in 1:num_cells(g)
        nodes = cell_nodes(g, i)
        @test nodes isa NTuple{4, Int}
        for n in nodes
            p = node_coordinates(g, n)
            @test abs(norm(p) - g.R) < 1e-10
        end
    end
end
