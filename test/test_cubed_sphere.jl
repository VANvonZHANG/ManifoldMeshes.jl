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

@testset "CubedSphereGrid geometry" begin
    g = CubedSphereGrid(n = 4)

    # Total area conservation
    total = sum(cell_volume(g, i) for i in 1:num_cells(g))
    @test total ≈ 4π * g.R^2 rtol = 1e-10

    # All centroids on sphere surface
    for i in 1:num_cells(g)
        c = cell_centroid(g, i)
        @test abs(norm(c) - g.R) < 1e-10
    end

    # All nodes on sphere surface
    for i in 1:num_nodes(g)
        p = node_coordinates(g, i)
        @test abs(norm(p) - g.R) < 1e-10
    end
end

@testset "CubedSphereGrid patch queries" begin
    g = CubedSphereGrid(n = 3)

    # Each face should have exactly n*n cells
    for face in 1:6
        cells_in_face = count(i -> cell_face(g, i) == face, 1:num_cells(g))
        @test cells_in_face == g.n * g.n
    end

    # cell_local_2d within bounds
    for i in 1:num_cells(g)
        li, lj = cell_local_2d(g, i)
        @test 1 <= li <= g.n
        @test 1 <= lj <= g.n
    end

    # Face 1 cells have local_2d covering full range
    face1_cells = [i for i in 1:num_cells(g) if cell_face(g, i) == 1]
    locals = [cell_local_2d(g, c) for c in face1_cells]
    @test length(unique(locals)) == g.n * g.n
end

@testset "CubedSphereGrid rotation" begin
    # 90-degree rotation around Z axis
    rot = SMatrix{3,3}(0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
    g = CubedSphereGrid(n = 2, rotation = rot)

    # A known point should be rotated
    p = node_coordinates(g, 1)
    # Original +Z face corner at (s=-1,t=-1) → (-1,-1,1) normalized
    # After rotation: (1, -1, 1) normalized — on +Z face but different position
    @test abs(norm(p) - 1.0) < 1e-10
end

@testset "CubedSphereGrid constructor validation" begin
    @test_throws ArgumentError CubedSphereGrid(n = 0)
    @test_throws ArgumentError CubedSphereGrid(n = 2, projection = :invalid)
    @test_throws ArgumentError CubedSphereGrid(n = 2, R = 0.0)
    @test_throws ArgumentError CubedSphereGrid(n = 2, R = -1.0)
end

@testset "CubedSphereGrid area conservation by projection" begin
    for proj in (:gnomonic, :equiangular)
        g = CubedSphereGrid(n = 4, projection = proj)
        total = sum(cell_volume(g, i) for i in 1:num_cells(g))
        @test total ≈ 4π * g.R^2 rtol = 1e-10
    end
end

@testset "CubedSphereGrid cell_nodes" begin
    g = CubedSphereGrid(n = 2)

    # cell_nodes returns NTuple{4,Int}
    nodes = cell_nodes(g, 1)
    @test nodes isa NTuple{4,Int}

    # All referenced nodes are on sphere surface
    for cell_id in 1:num_cells(g)
        for node_id in cell_nodes(g, cell_id)
            p = node_coordinates(g, node_id)
            @test abs(norm(p) - g.R) < 1e-10
        end
    end
end

@testset "CubedSphereGrid bounds checking" begin
    g = CubedSphereGrid(n = 2)

    @test_throws BoundsError node_coordinates(g, 0)
    @test_throws BoundsError node_coordinates(g, num_nodes(g) + 1)
    @test_throws BoundsError cell_volume(g, 0)
    @test_throws BoundsError cell_volume(g, num_cells(g) + 1)
    @test_throws BoundsError cell_centroid(g, 0)
    @test_throws BoundsError cell_centroid(g, num_cells(g) + 1)
    @test_throws BoundsError cell_nodes(g, 0)
    @test_throws BoundsError cell_nodes(g, num_cells(g) + 1)
end

@testset "CubedSphereGrid unimplemented stubs" begin
    g = CubedSphereGrid(n = 2)

    @test_throws ErrorException cell_cells(g, 1)
    @test_throws ErrorException node_cells(g, 1)
    @test_throws ErrorException cell_edges(g, 1)
    @test_throws ErrorException edge_length(g, 1)
    @test_throws ErrorException edge_midpoint(g, 1)
    @test_throws ErrorException edge_outward_normal(g, 1, 1)
end

@testset "CubedSphereGrid rotation effect" begin
    g0 = CubedSphereGrid(n = 2)
    rot = SMatrix{3,3}(0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
    g = CubedSphereGrid(n = 2, rotation = rot)

    # Rotation should change node positions
    p0 = node_coordinates(g0, 1)
    p = node_coordinates(g, 1)
    @test p != p0
    @test norm(p) ≈ norm(p0)
end

@testset "CubedSphereGrid radius scaling" begin
    R = 2.5
    g = CubedSphereGrid(n = 4, R = R)
    total = sum(cell_volume(g, i) for i in 1:num_cells(g))
    @test total ≈ 4π * R^2 rtol = 1e-10
end

@testset "CubedSphereGrid boundary" begin
    g = CubedSphereGrid(n = 2)
    @test isempty(boundary_nodes(g, :default))
    @test isempty(boundary_edges(g, :default))
end

@testset "CubedSphereGrid num_nodes and num_edges" begin
    g = CubedSphereGrid(n = 3)
    @test num_nodes(g) == 6 * (3 + 1)^2
    @test num_edges(g) == 12 * 3 * (3 + 1)
end
