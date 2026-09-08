using ManifoldMeshes
using Test

# A cube as a quad mesh: 8 corners, 6 faces, 12 edges (Euler: 8 + 6 - 12 = 2).
# Corners are 1-based ids; faces list corners in cyclic boundary order.
# The split variant triangulates the +x face (rows 1-2), padding with -1:
# 7 faces, 13 edges (Euler: 8 + 7 - 13 = 2).
const CUBE_NODES = [(1, 1, 1), (1, 1, -1), (1, -1, -1), (1, -1, 1),
    (-1, 1, 1), (-1, 1, -1), (-1, -1, -1), (-1, -1, 1)]
const CUBE_FACES = [1 4 3 2;   # +x
                    5 6 7 8;   # -x
                    1 2 6 5;   # +y
                    4 3 7 8;   # -y
                    1 4 8 5;   # +z
                    2 3 7 6]   # -z
const CUBE_FACES_SPLIT = [1 4 3 -1;  # +x lower triangle
                          1 3 2 -1;  # +x upper triangle
                          5 6 7 8;
                          1 2 6 5;
                          4 3 7 8;
                          1 4 8 5;
                          2 3 7 6]
const KS_SPLIT = [3, 3, 4, 4, 4, 4, 4]

@testset "connectivity derivation" begin
    topo = ManifoldMeshes._derive_mesh_topology(CUBE_FACES, fill(4, 6))
    @test topo.n_nodes == 8
    @test topo.n_edges == 12
    @test topo.n_nodes + 6 - topo.n_edges == 2      # Euler characteristic (sphere)
    # closed mesh: every edge has exactly 2 adjacent cells
    @test all(e -> ManifoldMeshes.n_neighbors(topo.edge_cells, e) == 2, 1:12)
    # every cell has 4 nonzero neighbors
    @test all(c -> all(>(0), topo.cell_cells[c]), 1:6)
    # cube corners touch 3 faces and 3 edges
    @test all(n -> length(topo.node_cells[n]) == 3, 1:8)
    @test all(n -> length(topo.node_edges[n]) == 3, 1:8)
    # edge_nodes rows are the 12 unique undirected pairs
    pairs = Set{Tuple{Int, Int}}()
    for e in 1:12
        v = topo.edge_nodes[e]
        @test length(v) == 2
        push!(pairs, (min(v[1], v[2]), max(v[1], v[2])))
    end
    @test length(pairs) == 12
    # cell_edges is in cyclic corner order: edge k connects corner k and k+1
    for c in 1:6
        ce = topo.cell_edges[c]
        @test length(ce) == 4
        for k in 1:4
            v = topo.edge_nodes[ce[k]]
            lo, hi = minmax(CUBE_FACES[c, k], CUBE_FACES[c, mod1(k + 1, 4)])
            @test (v[1], v[2]) == (lo, hi)
        end
    end
end

@testset "connectivity derivation with mixed cells" begin
    topo = ManifoldMeshes._derive_mesh_topology(CUBE_FACES_SPLIT, KS_SPLIT)
    @test topo.n_nodes == 8
    @test topo.n_edges == 13          # Euler: V + F - E == 2 with F == 7
    @test topo.n_nodes + 7 - topo.n_edges == 2
    # the new diagonal edge (1, 3) exists exactly once
    diag = [e for e in 1:13
            if Tuple(sort(collect(topo.edge_nodes[e]))) == (1, 3)]
    @test length(diag) == 1
    # it is shared by exactly the two triangles (cells 1 and 2)
    @test sort(collect(topo.edge_cells[diag[1]])) == [1, 2]
    # variable arity rows: triangles have 3 edges/cells, quads 4
    @test length(topo.cell_edges[1]) == 3
    @test length(topo.cell_edges[3]) == 4
    # vertex 1 now belongs to 4 cells (2 triangles + +y + +z)
    @test length(topo.node_cells[1]) == 4

    # positive fill values must not inflate the node maps
    bigfill = copy(CUBE_FACES_SPLIT)
    bigfill[bigfill .== -1] .= 999_999
    tf = ManifoldMeshes._derive_mesh_topology(bigfill, KS_SPLIT)
    @test tf.n_nodes == 8
end

@testset "connectivity derivation with boundary" begin
    # open mesh: drop the +x face; all 12 edges survive (each still referenced
    # by its other face), but the 4 edges of +x now have only 1 cell
    open_topo = ManifoldMeshes._derive_mesh_topology(CUBE_FACES[2:end, :], fill(4, 5))
    @test open_topo.n_edges == 12
    @test count(e -> ManifoldMeshes.n_neighbors(open_topo.edge_cells, e) == 1, 1:12) == 4
    # the four neighbors of +x (+y, -y, +z, -z = sliced rows 2-5) each lose
    # exactly one neighbor -> one 0 sentinel; -x (sliced row 1) is untouched
    for r in 2:5
        @test count(==(0), open_topo.cell_cells[r]) == 1
    end
    @test count(==(0), open_topo.cell_cells[1]) == 0
end
