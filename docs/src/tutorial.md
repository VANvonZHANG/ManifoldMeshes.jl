# Tutorial

This five-stage tutorial walks you from creating your first grid to producing publication-quality visualizations.
Every grid type in ManifoldMeshes.jl shares the same interface, so the skills you learn here transfer directly from `LatLonGrid` to `CubedSphereGrid`, `HEALPixGrid`, and `ReducedGaussianGrid`.

```@raw html
 <div style="border-left: 4px solid #4a90d9; padding: 0.5em 1em; margin: 1em 0; background: #f0f6ff; border-radius: 4px;">
<strong>Prerequisites.</strong>
Install ManifoldMeshes.jl and add a Makie backend for the visualization stage:
<pre><code class="language-julia">using Pkg
Pkg.add("ManifoldMeshes")
Pkg.add("CairoMakie")   # or GLMakie for interactive 3D
</code></pre>
</div>
```

## Stage 1: Create Your First Grid

The simplest grid is a latitude-longitude grid (`LatLonGrid`).
You specify the latitude and longitude cell *edges* as vectors of Float64.
The edges must start at ``-90^\circ`` / ``0^\circ`` and end at ``90^\circ`` / ``360^\circ``.

```julia
using ManifoldMeshes

# 5 latitude bands × 10 longitude bands = 50 cells
grid = LatLonGrid(
    lat_edges = collect(range(-90.0, 90.0; length = 6)),   # 6 edges → 5 bands
    lon_edges = collect(range(0.0, 360.0; length = 11)),   # 11 edges → 10 bands
)

println("Cells: ", num_cells(grid))   # 50
println("Nodes: ", num_nodes(grid))   # 66
println("Edges: ", num_edges(grid))   # 105
```

```
Cells: 50
Nodes: 66
Edges: 105
```

**What happened?** The constructor pre-computed every node position, cell volume, and cell centroid and cached them internally.
All subsequent queries are O(1) lookups.

## Stage 2: Query Geometry

Each cell has a *volume* (solid angle on the unit sphere, in steradians) and a *centroid* (the Fréchet mean of its vertices on S², returned as a 3D Cartesian `SVector`).

```julia
# Volume of cell 1 (a polar cell near the south pole)
v = cell_volume(grid, 1)
println("Volume of cell 1: ", round(v; digits = 6))

# Centroid — a point on the unit sphere
c = cell_centroid(grid, 1)
println("Centroid of cell 1: ", round.(c; digits = 4))
```

```
Volume of cell 1: 0.109662
Centroid of cell 1: [-0.309, 0.0, -0.9511]
```

A fundamental sanity check: the sum of all cell volumes on the unit sphere must equal ``4\pi``.

```julia
total = sum(cell_volume(grid, i) for i in 1:num_cells(grid))
println("Total solid angle: ", round(total; digits = 6))
println("4π =              ", round(4π; digits = 6))
```

```
Total solid angle: 12.566371
4π =              12.566371
```

This conservation law holds for *every* grid type in this library.

## Stage 3: Explore Topology

Mesh topology tells you who is next to whom.
All topology functions return integer IDs that you can feed back into any geometry query.

```julia
# Which nodes define cell 25?
nodes_25 = cell_nodes(grid, 25)
println("Nodes of cell 25: ", nodes_25)

# Which cells neighbor cell 25?  (0 = boundary sentinel — no neighbor)
neighbors_25 = cell_cells(grid, 25)
println("Neighbors of cell 25: ", neighbors_25)

# Which cells share node 7?
cells_around_7 = node_cells(grid, 7)
println("Cells sharing node 7: ", cells_around_7)
```

```
Nodes of cell 25: (30, 31, 41, 42)
Neighbors of cell 25: (24, 26, 15, 35)
Cells sharing node 7: [6, 7, 16, 17]
```

### Walking a Neighbor Chain

Because `cell_cells` returns the neighbors of any cell, you can walk across the sphere step by step.
Here we start at cell 1 (near the south pole) and move north through the first neighbor that is not our previous cell:

```julia
path = [1]
prev = 0
for _ in 1:5
    nbrs = cell_cells(grid, path[end])
    # Pick the first non-zero neighbor that is not where we came from
    for n in nbrs
        if n != 0 && n != prev
            prev = path[end]
            push!(path, n)
            break
        end
    end
end
println("Walk path: ", path)
```

```
Walk path: [1, 2, 3, 4, 5, 6]
```

Edge connectivity is also available via `cell_edges(grid, i)`, which returns the edge IDs surrounding a cell.

```julia
edges_25 = cell_edges(grid, 25)
println("Edges of cell 25: ", edges_25)
```

```
Edges of cell 25: (54, 59, 105, 100)
```

## Stage 4: Try Another Grid Type

The same interface works for every grid.
Switching to a `HEALPixGrid` requires only a different constructor — all query functions are identical.

```julia
hp = HEALPixGrid(nside = 4)    # 12 × 4² = 192 cells

println("HEALPix cells: ", num_cells(hp))    # 192
println("Volume sum:    ", round(
    sum(cell_volume(hp, i) for i in 1:num_cells(hp)); digits = 6))
println("Cell 1 nodes:  ", cell_nodes(hp, 1))
```

```
HEALPix cells: 192
Volume sum:    12.566371
Cell 1 nodes:  (1, 2, 5, 4)
```

### Cubed-Sphere: Patch Topology

The `CubedSphereGrid` maps six cube faces onto the sphere.
It provides two extra functions — `cell_face` and `cell_local_2d` — that identify which face a cell lives on and its in-face (i, j) indices.

```julia
cs = CubedSphereGrid(n = 4)    # 6 × 4² = 96 cells

println("Cubed-sphere cells: ", num_cells(cs))   # 96
println("Face of cell 1:     ", cell_face(cs, 1))
println("Local (i,j) of 1:   ", cell_local_2d(cs, 1))
```

```
Cubed-sphere cells: 96
Face of cell 1:     1
Local (i,j) of 1:   (1, 1)
```

Cross-face neighbors are seamless — `cell_cells` returns the correct global cell ID even when cells sit on the boundary between two cube faces.

```julia
# Cell on the edge of face 1 — its neighbor belongs to another face
edge_cell = 4   # last column of first row on face 1
nbrs = cell_cells(cs, edge_cell)
println("Cell $edge_cell neighbors: ", nbrs)
for n in nbrs
    n == 0 && continue
    println("  neighbor $n is on face ", cell_face(cs, n))
end
```

```
Cell 4 neighbors: (3, 5, 13, 68)
  neighbor 3 is on face 1
  neighbor 5 is on face 1
  neighbor 13 is on face 1
  neighbor 68 is on face 2
```

### Quick Comparison

```julia
using ManifoldMeshes

grids = [
    "LatLon 5×10"    => LatLonGrid(lat_edges = collect(range(-90.0, 90.0; length = 6)),
                                    lon_edges = collect(range(0.0, 360.0; length = 11))),
    "HEALPix Nside=4" => HEALPixGrid(nside = 4),
    "CubedSphere n=4"  => CubedSphereGrid(n = 4),
    "ReducedGauss T42" => ReducedGaussianGrid(nlat = 42),
]

for (label, g) in grids
    total = sum(cell_volume(g, i) for i in 1:num_cells(g))
    println(rpad(label, 20), " cells=", rpad(num_cells(g), 6),
            " total vol=", round(total; digits = 4))
end
```

```
LatLon 5×10          cells=50     total vol=12.5664
HEALPix Nside=4      cells=192    total vol=12.5664
CubedSphere n=4      cells=96     total vol=12.5664
ReducedGauss T42     cells=8192   total vol=12.5664
```

Every grid conserves ``4\pi``.

## Stage 5: Visualize

ManifoldMeshes.jl uses Makie for 3D rendering.
Load a Makie backend *before* calling any plotting function.

### Wireframe Plot

```julia
using CairoMakie   # or: using GLMakie

grid = HEALPixGrid(nside = 4)
fig = plot_mesh(grid; show_nodes = true)
save("healpix_wireframe.png", fig)
```

`plot_mesh` renders cell edges as great-circle arcs on a semi-transparent sphere.
Useful keyword arguments:

| Keyword | Default | Effect |
|---------|---------|--------|
| `show_edges` | `true` | Draw cell boundaries |
| `show_nodes` | `false` | Draw node markers |
| `show_cell_ids` | `false` | Label each cell with its index |
| `show_node_ids` | `false` | Label each node with its index |
| `n_arc_points` | `50` | Points per geodesic arc (smoothness) |
| `figsize` | `(800, 600)` | Figure size in pixels |

### Filled Plot

`plot_mesh_filled` colors each cell as a filled polygon, which is useful for displaying scalar fields.

```julia
# Color cells by their face ID (cubed-sphere only)
cs = CubedSphereGrid(n = 4)
face_colors = i -> cell_face(cs, i)  # returns 1–6, auto-mapped to colors
fig = plot_mesh_filled(cs; color_by = face_colors)
save("cubed_sphere_faces.png", fig)
```

The `color_by` keyword accepts any function `Int -> Color` (or `nothing` for the default light-blue fill).
Combine with `show_edges = true` to overlay cell boundaries.

### Interactive Exploration

If you use `GLMakie` instead of `CairoMakie`, the plot opens in an interactive window where you can rotate, zoom, and pan the sphere freely — no fixed camera angles to configure.

```julia
using GLMakie

grid = ReducedGaussianGrid(nlat = 42)
fig = plot_mesh_filled(grid; show_edges = true)
display(fig)   # opens interactive window
```

### Summary

You have learned to:

1. **Create** grids from latitude/longitude edges, HEALPix Nside, cubed-sphere resolution, or Gaussian latitude count.
2. **Query geometry** — cell volumes, centroids, node coordinates — and verify the ``4\pi`` conservation law.
3. **Traverse topology** — neighbor chains, incident cells, edge connectivity.
4. **Switch grid types** without changing any downstream code, thanks to the unified `AbstractManifoldMesh` interface.
5. **Visualize** wireframes and filled plots with Makie.

### Further Exploration

For a complete side-by-side comparison of all four grid types, run the
example script included in the repository:

```bash
julia --project=. examples/grids_visualization.jl
```

This produces a publication-quality figure showing coarse and fine
resolutions of each grid type, useful for choosing the right grid for
your application.

The [API Reference](api.md) documents every function in detail.
For the mathematical background behind geodesic cell boundaries and spherical geometry, see the [Theory](theory.md) page.
