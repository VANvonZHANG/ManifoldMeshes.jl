# Visualization

ManifoldMeshes.jl provides built-in visualization for sphere meshes via
[Makie.jl](https://makie.juliaplots.org/).
Two levels of API are available:

- **High-level plotting** -- [`plot_mesh`](@ref) and [`plot_mesh_filled`](@ref) produce
  complete figures with a semi-transparent sphere background.
- **Low-level data extraction** -- [`node_points`](@ref), [`edge_segments`](@ref),
  [`cell_polygons`](@ref), and [`slerp`](@ref) give you the raw geometry for custom
  rendering with any Makie backend.

## Setup

Load a Makie backend **before** calling any plot function.
For static images (documents, CI):

```julia
using CairoMakie
using ManifoldMeshes
```

For interactive 3D rotation and zoom:

```julia
using GLMakie
using ManifoldMeshes
```

!!! note
    The plotting functions are defined in an extension module that requires CairoMakie.
    If you want interactive rendering, switch the extension to load GLMakie instead,
    or extract the raw data (see [Data Extraction](@ref)) and build your own scene.

## Basic Wireframe

[`plot_mesh`](@ref) draws the mesh edges as great-circle arcs on a semi-transparent
sphere.  By default, edges are drawn but nodes are hidden.

```julia
using CairoMakie
using ManifoldMeshes

grid = LatLonGrid(
    lat_edges = collect(-90.0:30.0:90.0),
    lon_edges = collect(0.0:45.0:360.0),
)

fig = plot_mesh(grid)
save("wireframe.png", fig)
```

### Optional keyword arguments

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `show_nodes` | `Bool` | `false` | Draw nodes as scatter markers |
| `show_edges` | `Bool` | `true` | Draw edge arcs |
| `show_cell_ids` | `Bool` | `false` | Label each cell with its index |
| `show_node_ids` | `Bool` | `false` | Label each node with its index |
| `n_arc_points` | `Int` | `50` | Points per great-circle arc segment |
| `figsize` | `Tuple{Int,Int}` | `(800, 600)` | Figure size in pixels |

Show nodes and cell IDs to inspect topology:

```julia
fig = plot_mesh(grid; show_nodes = true, show_cell_ids = true)
```

## Filled Cells

[`plot_mesh_filled`](@ref) renders each cell as a filled spherical polygon.
The optional `color_by` keyword lets you color cells by index, data field, or any
criterion.

!!! note "Implementation detail"
    `plot_mesh_filled` uses batched rendering internally: a single
    `mesh!(verts, faces)` call with per-vertex coloring, rather than
    drawing each cell as a separate polygon. This improves performance
    significantly for large grids (hundreds to thousands of cells).

```julia
fig = plot_mesh_filled(grid)
save("filled.png", fig)
```

### Coloring by cell index

`color_by` accepts a callable that takes a cell index (`Int`) and returns a
Makie-compatible color:

```julia
using CairoMakie

# Color by latitude band: cycle through a palette
colors = Makie.wong_colors()
fig = plot_mesh_filled(grid; color_by = i -> colors[mod1(i, length(colors))])
```

### Optional keyword arguments

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `show_edges` | `Bool` | `true` | Overlay wireframe edges |
| `color_by` | `Function` or `nothing` | `nothing` | Cell color mapper; `nothing` = lightblue |
| `n_arc_points` | `Int` | `50` | Points per great-circle arc segment |
| `figsize` | `Tuple{Int,Int}` | `(800, 600)` | Figure size in pixels |

## Data Extraction

If you need full control over the rendering, use the low-level functions to obtain
discretized geometry as plain arrays of `Point3f`.

### [`node_points`](@ref)

```julia
pts = node_points(grid)  # Vector{Point3f}, indexed by node ID
```

Returns all node positions as `GeometryBasics.Point3f`, one per node, in node-ID order.

### [`edge_segments`](@ref)

```julia
segs = edge_segments(grid; n_arc_points = 20)
```

Returns a `Vector{Vector{Point3f}}`, one entry per edge.  Each inner vector is a
discretized great-circle arc connecting the two endpoint nodes.  Increase
`n_arc_points` for smoother arcs.

```julia
# Custom rendering with GLMakie
using GLMakie

fig = Figure()
ax = LScene(fig[1, 1]; show_axis = false)

segs = edge_segments(grid; n_arc_points = 30)
for seg in segs
    length(seg) > 1 && lines!(ax, seg; color = :steelblue, linewidth = 1)
end

display(fig)
```

### [`cell_polygons`](@ref)

```julia
polys = cell_polygons(grid; n_arc_points = 20)
```

Returns a `Vector{Vector{Point3f}}`, one entry per cell.  Each polygon traces the
closed cell boundary via great-circle arcs along the four edges (south, east,
north reversed, west reversed).

```julia
polys = cell_polygons(grid; n_arc_points = 30)

fig = Figure()
ax = LScene(fig[1, 1]; show_axis = false)

for (i, poly) in enumerate(polys)
    length(poly) >= 3 || continue
    center = poly[1]
    for k in 2:(length(poly) - 1)
        mesh!(ax, [center, poly[k], poly[k + 1]];
            color = :lightblue, transparency = true)
    end
end

display(fig)
```

### [`cell_triangles`](@ref)

```julia
verts, faces = cell_triangles(grid)
```

Returns a shared-vertex triangle mesh suitable for `GeometryBasics.Mesh`.
Each cell is decomposed into K triangles fanning from the sphere-projected
 centroid to consecutive boundary-node pairs, where K is the number of nodes
bounding that cell.

```julia
verts, faces = cell_triangles(grid)

# Custom rendering with per-cell coloring
using GLMakie, GeometryBasics

fig = Figure()
ax = LScene(fig[1, 1]; show_axis = false)

mesh!(ax, GeometryBasics.Mesh(verts, faces);
    color = :lightblue, shading = NoShading)

display(fig)
```

This is the primitive used internally by [`plot_mesh_filled`](@ref) for
batched filled rendering. You can call it directly to build custom
visualizations with per-vertex or per-face color arrays.

## Spherical Interpolation

[`slerp`](@ref) (spherical linear interpolation) computes evenly spaced points along
the great-circle arc between two unit-sphere points.

```julia
using ManifoldMeshes
using StaticArrays

p1 = SVector(1.0, 0.0, 0.0)
p2 = SVector(0.0, 1.0, 0.0)

# 10 points along the quarter-circle from p1 to p2
arc = slerp(p1, p2, 10)
```

This is the primitive used internally by [`edge_segments`](@ref) and
[`cell_polygons`](@ref) to discretize edges.  You can also call it directly to draw
custom arcs, geodesic paths, or interpolated trajectories on the sphere.

**Special cases:**
- Coincident endpoints (`theta < 1e-10`): returns a single-point vector.
- Antipodal endpoints (`theta ~ pi`): picks a stable intermediate great-circle plane
  and returns the half-great-circle path.

## Saving Figures

Use the standard Makie `save` function to write figures to disk:

```julia
fig = plot_mesh(grid)
save("mesh.png", fig)       # PNG raster
save("mesh.svg", fig)       # SVG vector (CairoMakie only)
```

Supported formats depend on the backend.  CairoMakie supports PNG, SVG, and PDF.
GLMakie supports PNG and JPEG.

## Example Script

The repository includes a standalone example script that produces a 4×2
comparison figure of all grid types at coarse and fine resolutions:

```bash
julia --project=. examples/grids_visualization.jl
```

This script demonstrates custom rendering with `edge_segments` and
[`LScene`](https://docs.makie.org/stable/api/#Makie.LScene) for full
control over layout and styling.
