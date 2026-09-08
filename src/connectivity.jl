# Generic mesh-topology derivation from a face-node table.
#
# Given a padded face-node table (1-based node ids, one row per cell, the
# first ks[c] entries of row c active and in cyclic boundary order), derives
# every adjacency CSR the mesh interface needs. Used by UnstructuredMesh; kept
# mesh-type-agnostic so future types and tests can reuse it. Edge IDs are
# assigned by first appearance while scanning cells in ascending order —
# deterministic across sessions.

"""
    _derive_mesh_topology(face_nodes::AbstractMatrix{Int}, ks::Vector{Int}) -> NamedTuple

Derive full polygon-mesh topology from a padded face-node table. `ks[c]` is the
number of active corners in row `c` (the leading non-fill entries); cells may
be triangles, quads, hexagons, ... (variable per cell).

Returns `(n_nodes, n_edges, cell_edges, edge_nodes, edge_cells, cell_cells,
node_cells, node_edges)`. All maps are `CSRMapping`s. `cell_edges[c][k]` is the
edge between corner k and corner k+1 of cell c (cyclic); `edge_cells` rows have
length 1 (mesh boundary) or 2 (interior); `cell_cells` uses 0 as the sentinel
for missing neighbors.
"""
function _derive_mesh_topology(face_nodes::AbstractMatrix{Int}, ks::Vector{Int})
    n_cells = size(face_nodes, 1)
    n_cells > 0 || throw(ArgumentError("face_nodes must have at least one row"))
    length(ks) == n_cells ||
        throw(ArgumentError("ks has $(length(ks)) entries for $n_cells cells"))
    n_nodes = maximum(face_nodes)

    # cell -> edges in cyclic corner order; edges keyed by undirected pair
    edge_ids = Dict{Tuple{Int, Int}, Int}()
    cell_edge_counts = ks
    _cell_edges, ceptrs = CSRMapping(n_cells, cell_edge_counts)
    for c in 1:n_cells
        for k in 1:ks[c]
            n1 = face_nodes[c, k]
            n2 = face_nodes[c, mod1(k + 1, ks[c])]
            e = get!(edge_ids, minmax(n1, n2)) do
                length(edge_ids) + 1
            end
            _cell_edges.values[ceptrs[c]] = e
            ceptrs[c] += 1
        end
    end
    n_edges = length(edge_ids)

    # edge -> nodes (ascending pair order)
    _edge_nodes = CSRMapping(n_edges, 2)
    for (key, e) in edge_ids
        _edge_nodes.values[_edge_nodes.offsets[e] - 1 + 1] = key[1]
        _edge_nodes.values[_edge_nodes.offsets[e] - 1 + 2] = key[2]
    end

    # edge -> cells (two-pass build; 1 or 2 cells per edge)
    edge_cell_counts = fill(0, n_edges)
    for c in 1:n_cells
        for k in 1:ks[c]
            e = _cell_edges.values[_cell_edges.offsets[c] - 1 + k]
            edge_cell_counts[e] += 1
        end
    end
    _edge_cells, ptrs = CSRMapping(n_edges, edge_cell_counts)
    for c in 1:n_cells
        for k in 1:ks[c]
            e = _cell_edges.values[_cell_edges.offsets[c] - 1 + k]
            _edge_cells.values[ptrs[e]] = c
            ptrs[e] += 1
        end
    end

    # cell -> cells via shared edges (0 sentinel at mesh boundary)
    _cell_cells, c2ptrs = CSRMapping(n_cells, ks)
    for c in 1:n_cells
        for k in 1:ks[c]
            e = _cell_edges.values[_cell_edges.offsets[c] - 1 + k]
            row = _edge_cells[e]
            other = length(row) == 2 ? (row[1] == c ? row[2] : row[1]) : 0
            _cell_cells.values[c2ptrs[c]] = other
            c2ptrs[c] += 1
        end
    end

    # node -> cells (two-pass)
    node_cell_counts = fill(0, n_nodes)
    for c in 1:n_cells, k in 1:ks[c]

        node_cell_counts[face_nodes[c, k]] += 1
    end
    _node_cells, ptrs2 = CSRMapping(n_nodes, node_cell_counts)
    for c in 1:n_cells, k in 1:ks[c]

        n = face_nodes[c, k]
        _node_cells.values[ptrs2[n]] = c
        ptrs2[n] += 1
    end

    # node -> edges (two-pass, from the edge->node table)
    node_edge_counts = fill(0, n_nodes)
    for e in 1:n_edges
        node_edge_counts[_edge_nodes.values[_edge_nodes.offsets[e] - 1 + 1]] += 1
        node_edge_counts[_edge_nodes.values[_edge_nodes.offsets[e] - 1 + 2]] += 1
    end
    _node_edges, ptrs3 = CSRMapping(n_nodes, node_edge_counts)
    for e in 1:n_edges
        n1 = _edge_nodes.values[_edge_nodes.offsets[e] - 1 + 1]
        n2 = _edge_nodes.values[_edge_nodes.offsets[e] - 1 + 2]
        _node_edges.values[ptrs3[n1]] = e
        ptrs3[n1] += 1
        _node_edges.values[ptrs3[n2]] = e
        ptrs3[n2] += 1
    end

    return (n_nodes = n_nodes, n_edges = n_edges, cell_edges = _cell_edges,
        edge_nodes = _edge_nodes, edge_cells = _edge_cells, cell_cells = _cell_cells,
        node_cells = _node_cells, node_edges = _node_edges)
end
