# HEALPixGrid

`HEALPixGrid` is a semi-structured grid on the unit sphere $S^2$ that partitions the sphere into cells of **exactly equal area**.
The name stands for Hierarchical Equal Area iso-Latitude Pixelization: the sphere is first divided into 12 base pixels, each of which is recursively subdivided into `nside` x `nside` quadrilateral cells, yielding a total of $12 \times \mathrm{nside}^2$ cells.
This grid is the standard pixelization in astrophysics (CMB maps, all-sky surveys) and provides a natural multi-resolution hierarchy.

```@raw html
<table style="margin: 1em 0;">
<tr><th style="text-align:left; padding-right:1.5em;">Trait</th><th style="text-align:left;">Value</th></tr>
<tr><td style="padding-right:1.5em;"><code>TopologyStyle</code></td><td><code>IsSemiGrid()</code></td></tr>
<tr><td style="padding-right:1.5em;"><code>CellTypeStyle</code></td><td><code>IsUniform{4}()</code></td></tr>
<tr><td style="padding-right:1.5em;"><code>PatchStyle</code></td><td><code>NoPatch()</code></td></tr>
</table>
```

## Construction

```julia
HEALPixGrid(; nside, ordering=:ring, rotation=I, R=1.0)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `nside` | `Int` | Resolution parameter. Must be >= 1. Total cells = $12 \times \mathrm{nside}^2$. |
| `ordering` | `Symbol` | Cell ordering scheme: `:ring` (default) or `:nested`. |
| `rotation` | `SMatrix{3,3,Float64}` | Rotation matrix applied to all node positions. Defaults to the identity. |
| `R` | `Float64` | Sphere radius (default `1.0`). Must be positive. |

In **ring** ordering, cells are numbered by iso-latitude rings from the north pole to the south pole.
In **nested** ordering, cells are numbered by base pixel (0--11), then by Morton (Z-order) index within each base pixel.
The nested scheme supports hierarchical multi-resolution operations.

```@example healpix
using ManifoldMeshes

# nside=4: 12 * 16 = 192 cells
grid = HEALPixGrid(; nside=4)

println("Cells: ", num_cells(grid))
println("Nodes: ", num_nodes(grid))
println("Edges: ", num_edges(grid))
```

A coarser grid with nested ordering:

```@example healpix
nested = HEALPixGrid(; nside=2, ordering=:nested)
println("Cells: ", num_cells(nested))
println("Ordering: ", nested.ordering)
```

A higher-resolution grid for realistic applications:

```@example healpix
fine = HEALPixGrid(; nside=16)
println("Cells: ", num_cells(fine))
```

## Cell Geometry

### Shape

Each cell is a spherical quadrilateral with four vertices returned in the order **SW, SE, NE, NW** (southwest, southeast, northeast, northwest).
The sphere is divided into three zones:

- **North polar cap** (rings 1 to `nside`): ring $r$ contains $4r$ cells, tapering from 4 cells at the pole to $4 \times \mathrm{nside}$ at the equatorial boundary.
- **Equatorial belt** (rings `nside+1` to `3*nside`): each ring contains exactly $4 \times \mathrm{nside}$ cells.
- **South polar cap** (rings `3*nside+1` to `4*nside-1`): symmetric with the north polar cap.

Cell vertices are computed at construction by normalizing the average of adjacent cell centers that surround each corner, with pole points used directly for polar-cap cells.
Shared corners between adjacent cells are **deduplicated** at construction time via coordinate rounding to ~$10^{-12}$ precision.

### Volume

All cells have **exactly equal area** by construction:

$$A_{\text{cell}} = \frac{4\pi R^2}{12 \times \mathrm{nside}^2}$$

Cell volumes are computed at construction using l'Huilier's formula for spherical triangle area (each quadrilateral is split into two triangles), then cached for O(1) lookup.

```@example healpix
# All cells have equal volume
v1 = cell_volume(grid, 1)
v100 = cell_volume(grid, 100)
println("Volume of cell   1: ", round(v1; digits=6))
println("Volume of cell 100: ", round(v100; digits=6))
println("Identical? ", v1 == v100)
```

### Centroid

Cell centroids are the Frechet mean (intrinsic mean) of the four vertices on $S^2$, computed via `Manifolds.mean`.
They are returned as `SVector{3, Float64}` in Cartesian coordinates.

```@example healpix
c = cell_centroid(grid, 1)
println("Centroid of cell 1: ", round.(c; digits=4))
```

### Area Conservation

The sum of all cell volumes on the unit sphere equals $4\pi$:

```@example healpix
total = sum(cell_volume(grid, i) for i in 1:num_cells(grid))
println("Sum of volumes: ", round(total; digits=6))
println("4*pi:           ", round(4pi; digits=6))
```

## Topology and Connectivity

HEALPixGrid has `IsSemiGrid` topology: cells are arranged in iso-latitude rings, but ring widths vary (polar caps vs.\ equatorial belt).
All connectivity is pre-computed and stored at construction time.

### Cell-to-Cell Neighbors

[`cell_cells`](@ref) returns a 4-tuple of neighbor cell IDs, one per edge.
Each non-degenerate edge is shared by exactly two cells.
When a cell edge has coincident endpoints (degenerate edge), the neighbor across that edge is `0`.

```@example healpix
# Neighbors of cell 1 (north polar cap)
neighbors = cell_cells(grid, 1)
println("Neighbors of cell 1: ", neighbors)

# Neighbors of a mid-latitude cell
mid = div(num_cells(grid), 2)
neighbors_mid = cell_cells(grid, mid)
println("Neighbors of cell $mid: ", neighbors_mid)
```

### Cell-to-Node Mapping

[`cell_nodes`](@ref) returns four node IDs as `(SW, SE, NE, NW)`.
Shared corners between adjacent cells have been merged to single node IDs at construction time.

```@example healpix
println("Nodes of cell 1: ", cell_nodes(grid, 1))
println("Nodes of cell $mid: ", cell_nodes(grid, mid))
```

### Node-to-Cell Incidence

[`node_cells`](@ref) returns a vector of cell IDs that share a given node.
Interior nodes typically belong to four cells.

```@example healpix
# How many cells share node 1?
nc = node_cells(grid, 1)
println("Cells around node 1: ", length(nc))
```

### Cell-to-Edge Mapping

[`cell_edges`](@ref) returns four edge IDs corresponding to the four sides of the quadrilateral.
Edges are derived from cell-node connectivity at construction: each unique unordered pair of adjacent nodes defines one edge.

### Boundary

Because HEALPixGrid covers the full sphere, [`boundary_nodes`](@ref) and [`boundary_edges`](@ref) always return empty vectors for any marker.

```@example healpix
println("Boundary nodes: ", boundary_nodes(grid, :default))
println("Boundary edges: ", boundary_edges(grid, :default))
```

## When to Use

### Strengths

- **Exact equal area.** Every cell has the same solid angle $4\pi / (12 \times \mathrm{nside}^2)$, making HEALPix ideal for Monte Carlo integration, histogram-based statistics, and any application where uniform sampling is required.
- **Hierarchical structure.** The nested ordering and power-of-2 `nside` values create a natural multi-resolution pyramid: coarsening from `nside` to `nside/2` is a simple pixel aggregation.
- **Iso-latitude rings.** Cells are organized along rings of constant latitude, enabling fast spherical-harmonic transforms via the HEALPix library's ring-based algorithms.
- **No polar singularity.** Unlike [LatLonGrid](latlon.md), cells near the poles have the same area and shape as equatorial cells.

### Weaknesses

- **Semi-structured indexing.** There is no simple Cartesian `(i, j)` index in ring ordering; cell-to-ring mapping requires accumulated ring counts.
- **Cell shapes vary.** Although areas are equal, cell shapes are not: polar-cap cells are more elongated than equatorial cells. This matters for finite-volume schemes where cell aspect ratio affects numerical accuracy.
- **Resolution constrained to powers of 2.** For full hierarchical support, `nside` should be a power of 2. Non-power-of-2 values work for ring ordering but break the nested hierarchy.

### Typical Use Cases

- Cosmic Microwave Background (CMB) map analysis and visualization
- All-sky astronomical surveys (Gaia, WMAP, Planck)
- Statistical analysis on the sphere requiring uniform sampling (pixel statistics, angular power spectra)
- Hierarchical spatial indexing for large-scale point catalogs
- Multi-resolution analysis on the sphere (wavelet transforms, needlet filtering)

## Gotchas

1. **Cell IDs depend on ordering.** Ring and nested ordering assign different cell IDs to the same physical region of the sphere. Code that stores or compares cell IDs must be consistent about which ordering it expects.

2. **`node_cells` is O(N).** The current implementation scans all cells to find those adjacent to a given node. For large grids (e.g.\ `nside=1024` gives 12,582,912 cells), repeated calls to [`node_cells`](@ref) can be expensive. Consider caching results if you need them frequently.

3. **Polar-cap cell shapes are elongated.** Although all cells have equal area, polar-cap cells (rings near the poles) have a higher aspect ratio than equatorial cells. This is inherent to the HEALPix design and is the trade-off for iso-latitude ring structure combined with equal area.

4. **No boundary.** Because the grid covers the full closed sphere, [`boundary_nodes`](@ref) and [`boundary_edges`](@ref) return empty vectors for any marker argument. There is no partial-sphere or regional mode.

5. **`edge_outward_normal` return type.** The function returns a `NamedTuple{(:base_point, :normal)}` (tangent-space semantics), not a plain vector. For degenerate zero-length edges, the normal is `zero(SVector{3,Float64})`.

6. **Corner deduplication uses coordinate rounding.** Shared nodes between adjacent cells are merged by rounding coordinates to $10^{-12}$ precision. In theory, two very close but distinct corners could be incorrectly merged, though this is extremely unlikely at any practical `nside` value.
