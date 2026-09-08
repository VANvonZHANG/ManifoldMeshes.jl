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

using LinearAlgebra   # normalize/norm in the R=3 analytic block
using ManifoldMeshes: LatLonGrid, cell_centroid, cell_volume, edge_cells,
                      edge_length, edge_midpoint, edge_outward_normal, edge_nodes

# Build an UnstructuredMesh carrying the same nodes/connectivity as a LatLonGrid.
function unstructured_from_grid(g)
    fn = Matrix{Int}(undef, num_cells(g), 4)
    for c in 1:num_cells(g)
        fn[c, :] .= collect(cell_nodes(g, c))
    end
    node_lon = Vector{Float64}(undef, num_nodes(g))
    node_lat = Vector{Float64}(undef, num_nodes(g))
    for n in 1:num_nodes(g)
        lat, lon = ManifoldMeshes._cartesian_to_latlon(node_coordinates(g, n))
        node_lon[n] = lon
        node_lat[n] = lat
    end
    return UnstructuredMesh(node_lon, node_lat, fn; R = g.R, start_index = 1)
end

@testset "geometry oracle vs LatLonGrid" begin
    # Full-span grid (LatLonGrid requires -90..90). Polar rows are included:
    # the fan triangulation from corner 1 equals LatLonGrid's A-C diagonal
    # split bit-for-bit, including the degenerate polar triangles.
    g = LatLonGrid(lat_edges = collect(range(-90.0, 90.0; length = 17)),
        lon_edges = collect(range(0.0, 360.0; length = 33)))
    m = unstructured_from_grid(g)
    @test num_cells(m) == num_cells(g)
    @test num_nodes(m) == num_nodes(g)
    # LatLonGrid stores vertical edges only per left column (its east seam
    # edge maps to the column-1 edge); the face-node table references the
    # lon=360 seam nodes, so the derived mesh carries nlat extra vertical
    # seam edges — g's edge-key set is a strict subset of m's.
    @test num_edges(m) == num_edges(g) + 16
    for c in 1:num_cells(g)
        @test collect(cell_nodes(m, c)) == collect(cell_nodes(g, c))
        @test cell_volume(m, c) ≈ cell_volume(g, c)                # same fan split
        @test isapprox(cell_centroid(m, c), cell_centroid(g, c))   # both R=1 Riemannian means
    end

    # Edge ids enumerate differently (scan order vs LatLon's h-then-v); match by node pair.
    edge_key(mesh, e) = Tuple(sort(collect(edge_nodes(mesh, e))))
    g_edge_map(gg) = Dict(edge_key(gg, e) => e for e in 1:num_edges(gg))
    gem = g_edge_map(g)
    for e in [1, 37, 200, num_edges(m)]
        ge = gem[edge_key(m, e)]
        @test edge_length(m, e) ≈ edge_length(g, ge)
        @test edge_nodes(m, e) isa NTuple{2, Int}
        @test isapprox(edge_midpoint(m, e), edge_midpoint(g, ge))
        cell = edge_cells(m, e)[1]
        gcell = edge_cells(g, ge)[1]
        @test cell == gcell
        nm = edge_outward_normal(m, e, cell)
        ng = edge_outward_normal(g, ge, gcell)
        @test isapprox(nm.base_point, ng.base_point)
        @test isapprox(nm.normal, ng.normal; atol = 1e-12)
    end

    # non-unit radius: unit-normalized math keeps edge geometry correct —
    # checked ANALYTICALLY on a cube mesh at R = 3 (LatLonGrid's own edge math
    # goes through the unit Sphere(2) manifold, so it is not a valid oracle
    # at R != 1)
    m3 = UnstructuredMesh(Sphere(2), cube_points(3.0), CUBE_FACES)
    @test m3.R ≈ 3.0
    e14 = findfirst(e -> Tuple(sort(collect(edge_nodes(m3, e)))) == (1, 4),
        1:num_edges(m3))
    @test edge_length(m3, e14) ≈ 3.0 * acos(1 / 3)
    mp = edge_midpoint(m3, e14)
    @test norm(mp) ≈ 3.0
    @test isapprox(mp, 3.0 * normalize(SVector(1.0, 0.0, 1.0)); atol = 1e-12)
    # the six cube faces are congruent: each covers 4π/6 steradians
    @test all(c -> isapprox(cell_volume(m3, c), 9 * 4π / 6; rtol = 1e-12), 1:6)
end

@testset "geometry on mixed cells" begin
    quad_cube = UnstructuredMesh(Sphere(2), cube_points(), CUBE_FACES)
    split_cube = UnstructuredMesh(Sphere(2), cube_points(), CUBE_FACES_SPLIT;
        fill_value = -1)
    # the two triangles of the split +x face sum to the quad's area (fan from
    # vertex 1 uses the SAME 1-3 diagonal the quad's A-C split uses)
    @test cell_volume(split_cube, 1) + cell_volume(split_cube, 2) ≈
          cell_volume(quad_cube, 1) rtol = 1e-12
end
