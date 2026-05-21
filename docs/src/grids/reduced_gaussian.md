# ReducedGaussianGrid

`ReducedGaussianGrid` is a semi-structured grid on the unit sphere $S^2$ that uses Gaussian latitude bands with progressively fewer longitude cells toward the poles.
The longitude reduction follows an octahedral pattern, producing near-equal cell areas without the polar singularity of a regular lat-lon grid.
This grid type is the standard discretization for spectral weather and climate models (IFS, GFS).

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
ReducedGaussianGrid(; nlat, R=1.0)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `nlat` | `Int` | Number of latitude cell bands. Must be at least `2`. |
| `R` | `Float64` | Sphere radius (default `1.0`). Must be positive. |

The grid places `nlat + 1` latitude circles from the south pole to the north pole.
The interior `nlat - 1` boundaries are placed at equal-area intervals: $\sin(\theta_j) = -1 + 2j/n$ for $j = 1, \ldots, n{-}1$.
On each latitude circle, the number of longitude nodes follows an octahedral rule that increases from the poles toward the equator, keeping cell areas approximately uniform.

The total number of cells depends on `nlat` via the octahedral longitude counts and is not a simple closed-form expression.

```@example rgg
using ManifoldMeshes

# 8 latitude bands with octahedral longitude reduction
grid = ReducedGaussianGrid(; nlat=8)

println("Cells: ", num_cells(grid))
println("Nodes: ", num_nodes(grid))
println("Edges: ", num_edges(grid))
```

A coarser grid for quick experiments:

```@example rgg
coarse = ReducedGaussianGrid(; nlat=4)
println("Cells: ", num_cells(coarse))
```

## Cell Geometry

### Shape

Each cell is a spherical quadrilateral with four vertices returned in the order **SW, SE, NE, NW** (southwest, southeast, northeast, northwest).
Cells at the first and last latitude bands (adjacent to the poles) may degenerate into triangles because the pole node is shared; these are still represented as quadrilaterals with repeated (coincident) node IDs.

The vertex positions are projected onto $S^2$ via:

$$\mathbf{x} = R\,(\sin\theta\cos\phi,\;\sin\theta\sin\phi,\;\cos\theta)$$

where $\theta$ is the colatitude and $\phi$ is the longitude.

### Volume

Cell volumes are computed at construction using l'Huilier's formula for spherical triangle area. Each quadrilateral is split into two triangles along a diagonal (A-C), and their areas are summed. All volumes are cached for O(1) lookup.

```@example rgg
# Volume of cell 1
v1 = cell_volume(grid, 1)
println("Volume of cell 1: ", round(v1; digits=6))

# Volume of a mid-latitude cell
v_mid = cell_volume(grid, div(num_cells(grid), 2))
println("Volume of mid cell: ", round(v_mid; digits=6))
```

### Centroid

Cell centroids are the Frechet mean (intrinsic mean) of the four vertices on $S^2$, computed via `Manifolds.mean`.
They are returned as `SVector{3, Float64}` in Cartesian coordinates.

```@example rgg
c = cell_centroid(grid, 1)
println("Centroid of cell 1: ", round.(c; digits=4))
```

### Area Conservation

The sum of all cell volumes on the unit sphere equals $4\pi$:

```@example rgg
total = sum(cell_volume(grid, i) for i in 1:num_cells(grid))
println("Sum of volumes: ", round(total; digits=6))
println("4*pi:           ", round(4pi; digits=6))
```

## Topology and Connectivity

ReducedGaussianGrid has `IsSemiGrid` topology: cells are arranged in latitude bands, but each band contains a different number of cells.
There is no simple Cartesian indexing; instead, all connectivity is pre-computed and stored at construction time.

### Cell-to-Cell Neighbors

[`cell_cells`](@ref) returns a 4-tuple of neighbor cell IDs, one per edge.
When a cell edge is degenerate (zero-length, which can occur at the poles from duplicated nodes), the neighbor across that edge is `0`.
Otherwise, each non-degenerate edge is shared by exactly two cells.

```@example rgg
# Check neighbors of a cell near the equator
mid = div(num_cells(grid), 2)
neighbors = cell_cells(grid, mid)
println("Neighbors of cell $mid: ", neighbors)
```

### Cell-to-Node Mapping

[`cell_nodes`](@ref) returns four node IDs as `(SW, SE, NE, NW)`.
For polar cells, some of these IDs may be repeated (coincident nodes at the pole).

```@example rgg
# Nodes of cell 1 (near south pole)
println("Nodes of cell 1: ", cell_nodes(grid, 1))

# Nodes of a mid-latitude cell
println("Nodes of mid cell: ", cell_nodes(grid, mid))
```

### Node-to-Cell Incidence

[`node_cells`](@ref) returns a vector of cell IDs that share a given node.
Interior nodes typically belong to four cells; pole nodes may be shared by many cells.

```@example rgg
# How many cells share node 1?
nc = node_cells(grid, 1)
println("Cells around node 1: ", length(nc))
```

### Cell-to-Edge Mapping

[`cell_edges`](@ref) returns four edge IDs corresponding to the four sides of the quadrilateral.
Edges are derived from cell-node connectivity: each unique unordered pair of adjacent nodes defines one edge.

### Boundary

Because ReducedGaussianGrid covers the full sphere, [`boundary_nodes`](@ref) and [`boundary_edges`](@ref) always return empty vectors for any marker.

```@example rgg
println("Boundary nodes: ", boundary_nodes(grid, :default))
println("Boundary edges: ", boundary_edges(grid, :default))
```

## When to Use

### Strengths

- **Near-equal cell areas.** The octahedral longitude reduction keeps cell sizes roughly uniform from equator to pole, avoiding the extreme area variation of [LatLonGrid](latlon.md).
- **No polar singularity.** Unlike regular lat-lon grids, longitude lines do not converge at the poles, so there is no CFL restriction or numerical conditioning problem near the poles.
- **Spectral model compatibility.** This is the standard grid for spectral transform methods (Gaussian latitudes are the natural quadrature points for Legendre transforms). Data on this grid maps directly to spherical harmonics.
- **Pre-computed connectivity.** All topology (edges, neighbors) is derived at construction and cached, giving O(1) lookups.

### Weaknesses

- **Semi-structured indexing.** Because each latitude band has a different number of cells, there is no simple `(i, j)` Cartesian index. Accessing cells by latitude band requires knowing the cumulative cell offset.
- **Approximate equal area.** Cell areas are near-equal but not exactly equal (unlike [HEALPixGrid](healpix.md), which guarantees exact equal area).
- **Less intuitive than LatLonGrid.** The octahedral longitude pattern and variable band widths make this grid harder to reason about manually than a simple lat-lon grid.

### Typical Use Cases

- Spectral weather and climate models (IFS, GFS, MPAS)
- Gaussian quadrature on the sphere for numerical integration
- Situations where near-uniform resolution is needed but the full equal-area guarantee of HEALPix is not required
- Interfacing with existing spectral model output or initial conditions

## Gotchas

1. **Polar cells may have repeated node IDs.** Cells adjacent to the poles can have two or more identical node IDs because the single pole node is reused. This means [`cell_nodes`](@ref) may return a tuple like `(a, b, c, c)` where the last two IDs are the same point. Topology queries handle this correctly (degenerate edges get neighbor ID `0`), but comparing node IDs directly will reveal the duplication.

2. **Zero-length edges at the poles.** Edges between duplicated pole nodes have zero length. [`edge_length`](@ref) returns `0.0` for these edges, [`edge_midpoint`](@ref) returns the pole coordinate, and [`edge_outward_normal`](@ref) returns `(base_point=p, normal=zero(SVector{3,Float64}))`.

3. **`cell_cells` returns `0` for degenerate edges.** When an edge has coincident endpoints (self-loop), the neighbor across that edge is reported as `0`. Always check for `0` before using a neighbor ID as an array index.

4. **No boundary.** Because the grid covers the full closed sphere, [`boundary_nodes`](@ref) and [`boundary_edges`](@ref) return empty vectors for any marker argument. There is no partial-sphere or regional mode.

5. **`node_cells` is O(N).** The current implementation scans all cells to find those adjacent to a given node. For large grids, repeated calls to [`node_cells`](@ref) can be expensive. Consider caching results if you need them frequently.

6. **Area is near-equal, not exact.** Unlike [HEALPixGrid](healpix.md), ReducedGaussianGrid does not guarantee exactly equal cell areas. The octahedral reduction provides a good approximation, but cell areas can vary by up to roughly a factor of 2 between the smallest and largest cells depending on `nlat`.
