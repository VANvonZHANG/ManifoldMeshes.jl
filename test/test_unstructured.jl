using ManifoldMeshes
using Manifolds
using StaticArrays
using Test

using ManifoldMeshes: AbstractManifoldMesh, CellTypeStyle, IsMesh, IsMixed,
                      IsUniform, NoPatch, PatchStyle, TopologyStyle, cell_nodes,
                      manifold, num_cells, num_edges, num_nodes, node_coordinates

const CUBE_NODES = [(1, 1, 1), (1, 1, -1), (1, -1, -1), (1, -1, 1),
    (-1, 1, 1), (-1, 1, -1), (-1, -1, -1), (-1, -1, 1)]
const CUBE_FACES = [1 4 3 2; 5 6 7 8; 1 2 6 5; 4 3 7 8; 1 4 8 5; 2 3 7 6]
const CUBE_FACES_SPLIT = [1 4 3 -1; 1 3 2 -1; 5 6 7 8; 1 2 6 5; 4 3 7 8;
                          1 4 8 5; 2 3 7 6]

cube_points(R = 1.0) = [R / sqrt(3) * SVector{3, Float64}(c...) for c in CUBE_NODES]

@testset "UnstructuredMesh construction" begin
    pts = cube_points(2.5)
    m = UnstructuredMesh(Sphere(2), pts, CUBE_FACES)
    @test m isa AbstractManifoldMesh
    @test num_nodes(m) == 8
    @test num_cells(m) == 6
    @test num_edges(m) == 12
    @test TopologyStyle(m) === IsMesh()
    @test CellTypeStyle(m) === IsMixed{4}()
    @test PatchStyle(m) === NoPatch()
    @test manifold(m) == Sphere(2)
    @test m.R ≈ 2.5                                # inferred from mean node norm
    @test node_coordinates(m, 3) == pts[3]          # stored verbatim
    @test cell_nodes(m, 1) == (1, 4, 3, 2)
    @test cell_nodes(m, 6) == (2, 3, 7, 6)
    @test length(cell_edges(m, 1)) == 4

    # MixedCellTopology Tuple equality (direct, both directions)
    @test ManifoldMeshes.MixedCellTopology((1, 2, 3, 0), 3) == (1, 2, 3)
    @test (1, 2, 3) == ManifoldMeshes.MixedCellTopology((1, 2, 3, 0), 3)
    @test ManifoldMeshes.MixedCellTopology((1, 2, 3, 0), 3) != (1, 2, 4)
    @test ManifoldMeshes.MixedCellTopology((1, 2, 3, 0), 3) != (1, 2, 3, 0)

    # lon/lat convenience constructor with a 0-based (UGRID-style) table
    m2 = UnstructuredMesh(fill(10.0, 8), fill(20.0, 8), CUBE_FACES .- 1;
        R = 1.0, start_index = 0)
    @test num_cells(m2) == 6
    @test cell_nodes(m2, 1) == (1, 4, 3, 2)
    @test all(n -> node_coordinates(m2, n)[3] ≈ sind(20.0), 1:8)  # z = sin(lat)*R
end

@testset "UnstructuredMesh mixed cells" begin
    pts = cube_points()
    m = UnstructuredMesh(Sphere(2), pts, CUBE_FACES_SPLIT; fill_value = -1)
    @test CellTypeStyle(m) === IsMixed{4}()
    @test num_cells(m) == 7
    @test num_nodes(m) == 8
    @test num_edges(m) == 13                  # Euler: V + F - 2
    @test cell_nodes(m, 1) == (1, 4, 3)
    @test cell_nodes(m, 3) == (5, 6, 7, 8)
    @test length(cell_edges(m, 1)) == 3
    @test length(node_cells(m, 1)) == 4       # 2 triangles + +y + +z

    # single hexagon cell (polar hexagon at lat 30)
    m6 = UnstructuredMesh(collect(0.0:60.0:300.0), fill(30.0, 6),
        Matrix{Int}(reshape(1:6, 1, 6)); start_index = 1)
    @test CellTypeStyle(m6) === IsMixed{6}()
    @test num_cells(m6) == 1
    @test num_edges(m6) == 6
    @test cell_nodes(m6, 1) == (1, 2, 3, 4, 5, 6)
end

@testset "UnstructuredMesh validation" begin
    pts = cube_points()
    @test_throws ArgumentError UnstructuredMesh(Sphere(2), pts, ones(Int, 2, 2))
    @test_throws ArgumentError UnstructuredMesh(Sphere(2), pts[1:5],
        [1 2 3 4; 2 3 4 6])                       # node 6 out of 1:5 range
    @test_throws ArgumentError UnstructuredMesh(Sphere(2), pts[1:5],
        [1 2 3 3; 2 3 4 5])                       # repeated node in a cell
    @test_throws ArgumentError UnstructuredMesh(Sphere(2), pts,
        [1 2 -1 4; 2 3 4 5])                      # truncated prefix (2 active corners)
    @test_throws ArgumentError UnstructuredMesh(Sphere(2), pts,
        [1 2 -1 -1; 2 3 4 5])                     # only 2 active corners
    @test_throws ArgumentError UnstructuredMesh(Euclidean(2), pts,
        CUBE_FACES)                               # v1 gate: Sphere only
end
