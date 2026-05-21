# LatLonGrid

`LatLonGrid` is a structured latitude-longitude grid on the unit sphere $S^2$.
It is the simplest grid type in ManifoldMeshes.jl: cells are defined by the intersection of latitude bands and longitude sectors, producing a regular array of spherical quadrilaterals with Cartesian `(i, j)` indexing.
All node positions, cell volumes, and cell centroids are pre-computed at construction time, making every subsequent query an O(1) lookup.

```@raw html
<table style="margin: 1em 0;">
<tr><th style="text-align:left; padding-right:1.5em;">Trait</th><th style="text-align:left;">Value</th></tr>
<tr><td style="padding-right:1.5em;"><code>TopologyStyle</code></td><td><code>IsGrid()</code></td></tr>
<tr><td style="padding-right:1.5em;"><code>CellTypeStyle</code></td><td><code>IsUniform{4}()</code></td></tr>
<tr><td style="padding-right:1.5em;"><code>PatchStyle</code></td><td><code>NoPatch()</code></td></tr>
</table>
```

## Construction

```julia
LatLonGrid(; lat_edges, lon_edges, R=1.0)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `lat_edges` | `Vector{Float64}` | Latitude cell edges in degrees. Must start at ``-90^\circ`` and end at ``90^\circ``, monotonically increasing. |
| `lon_edges` | `Vector{Float64}` | Longitude cell edges in degrees. Must start at ``0^\circ`` and end at ``360^\circ``, monotonically increasing. |
| `R` | `Float64` | Sphere radius (default `1.0`). Must be positive. |

The number of cells is `(length(lat_edges) - 1) * (length(lon_edges) - 1)`.
The constructor copies the edge arrays internally to prevent external mutation from corrupting cached geometry.

```@example latlon
using ManifoldMeshes

# 4 latitude bands x 8 longitude bands = 32 cells
grid = LatLonGrid(
    lat_edges = collect(-90.0:45.0:90.0),   # 5 edges -> 4 bands
    lon_edges = collect(0.0:45.0:360.0),    # 9 edges -> 8 bands
)

println("Cells: ", num_cells(grid))   # 32
println("Nodes: ", num_nodes(grid))   # 45
println("Edges: ", num_edges(grid))   # 72
```

A convenience for creating regular grids with uniform spacing:

```@example latlon
# 1-degree resolution: 180 lat bands x 360 lon bands = 64 800 cells
fine = LatLonGrid(
    lat_edges = collect(-90.0:1.0:90.0),
    lon_edges = collect(0.0:1.0:360.0),
)

println("Cells: ", num_cells(fine))
```

## Cell Geometry

### Shape

Each cell is a spherical quadrilateral whose four vertices are the $(\theta, \phi)$ grid points projected onto $S^2$:

$$\mathbf{x} = R\,(\cos\theta\cos\phi,\;\cos\theta\sin\phi,\;\sin\theta)$$

Cell nodes are returned in the order **SW, SE, NE, NW** (southwest, southeast, northeast, northwest).

### Volume

Cell volumes (solid angles on the unit sphere, scaled by $R^2$) are computed at construction using l'Huilier's formula. Each quadrilateral is split into two spherical triangles along a diagonal. If the primary diagonal degenerates (coincident endpoints, which happens at the poles), the code falls back to the alternate diagonal. If both diagonals degenerate, a lune formula is used.

```@example latlon
# Volume of cell 1 (near the south pole)
v = cell_volume(grid, 1)
println("Volume of cell 1: ", round(v; digits=6))

# Volume of a mid-latitude cell (cell 9)
v_mid = cell_volume(grid, 9)
println("Volume of cell 9: ", round(v_mid; digits=6))
```

### Centroid

Cell centroids are the Fréchet mean (intrinsic mean) of the four vertices on $S^2$, computed via `Manifolds.mean`. They are returned as `SVector{3, Float64}` in Cartesian coordinates.

```@example latlon
c = cell_centroid(grid, 9)
println("Centroid of cell 9: ", round.(c; digits=4))
```

### Area Conservation

A fundamental sanity check: the sum of all cell volumes on the unit sphere equals $4\pi$.

```@example latlon
total = sum(cell_volume(grid, i) for i in 1:num_cells(grid))
println("Sum of volumes: ", round(total; digits=6))
println("4*pi:           ", round(4pi; digits=6))
```

## Topology and Connectivity

LatLonGrid uses Cartesian `(ilat, ilon)` indexing with a linear cell ID computed as:

```
cell_id = (ilat - 1) * nlon + ilon
```

### Cell-to-Cell Neighbors

[`cell_cells`](@ref) returns a 4-tuple `(south, north, west, east)` of neighbor cell IDs.
Longitude is periodic: the cell at `ilon = nlon` wraps to `ilon = 1` in the east direction.
At the north and south boundaries (first and last latitude bands), the missing neighbor is reported as `0`.

```@example latlon
# Cell at the equator, first longitude band
eq_cell = 9  # ilat=2, ilon=1 in this grid
south, north, west, east = cell_cells(grid, eq_cell)
println("Neighbors of cell $eq_cell: S=$south, N=$north, W=$west, E=$east")

# South-pole band cell (first lat band, ilat=1): south neighbor = 0
polar_cell = 1
south, north, west, east = cell_cells(grid, polar_cell)
println("Neighbors of cell $polar_cell: S=$south, N=$north, W=$west, E=$east")
```

### Cell-to-Node Mapping

[`cell_nodes`](@ref) returns the four corner node IDs as `(SW, SE, NE, NW)`.
The SE and NE indices wrap periodically when `ilon = nlon`.

```@example latlon
nodes = cell_nodes(grid, 1)
println("Nodes of cell 1: ", nodes)
```

### Node-to-Cell Incidence

[`node_cells`](@ref) returns a vector of cell IDs adjacent to a given node.
Interior nodes have four adjacent cells; edge and polar nodes have fewer.

```@example latlon
# South pole node (node 1)
nc = node_cells(grid, 1)
println("Cells around south pole: ", length(nc))
```

### Cell-to-Edge Mapping

[`cell_edges`](@ref) returns four edge IDs as `(south, north, west, east)`.
Edges are numbered in two blocks: horizontal edges first, then vertical edges.

### Boundary

Because LatLonGrid covers the full sphere, [`boundary_nodes`](@ref) and [`boundary_edges`](@ref) always return empty vectors for any marker.

```@example latlon
println("Boundary nodes: ", boundary_nodes(grid, :default))
println("Boundary edges: ", boundary_edges(grid, :default))
```

## When to Use

### Strengths

- **Simple structure.** Cartesian `(i, j)` indexing makes it easy to integrate with existing lat-lon datasets, NetCDF files, and legacy atmospheric models.
- **O(1) lookups.** All geometry is pre-computed; neighbor queries are pure arithmetic with no pointer chasing.
- **Periodic wrapping.** East-west periodicity is built in, so no special boundary handling is needed in the longitude direction.
- **Conservation.** The sum of all cell volumes equals $4\pi R^2$ to machine precision.

### Weaknesses

- **Polar singularity.** All longitude lines converge at the poles, causing cells to degenerate into extremely narrow wedges. This leads to the CFL restriction in explicit time-stepping and poor numerical conditioning near the poles.
- **Non-uniform cell area.** Cells near the equator are much larger than cells near the poles. A $1^\circ \times 1^\circ$ cell at the equator has roughly 57 times the area of one at $89^\circ$ latitude.
- **No quasi-uniform resolution.** If approximately equal cell areas are needed, consider [`CubedSphereGrid`](@ref), [`HEALPixGrid`](@ref), or [`ReducedGaussianGrid`](@ref) instead.

### Typical Use Cases

- Interfacing with legacy lat-lon data products (reanalysis, satellite level-3 binned data)
- Testing and prototyping (simplest grid to construct and reason about)
- Regional grid studies where polar degeneration is outside the domain of interest
- Reference solutions where the structured indexing simplifies analytic comparisons

## Gotchas

1. **Longitude 0 and longitude 360 are the same physical point** but have different linear node IDs. Node `(ilat, 1)` and node `(ilat, nlon+1)` coincide on the sphere. Topology functions handle this automatically via periodic wrapping, but if you compare node coordinates directly they will be equal while their IDs differ.

2. **Polar cells degenerate.** Cells adjacent to the poles have very small area, and edges along the pole (horizontal edges at `ilat=1` and `ilat=nlat+1`) collapse to zero length. [`edge_length`](@ref) returns `0.0` for these edges, and [`edge_midpoint`](@ref) returns one of the endpoints.

3. **`cell_cells` returns `0` at the poles.** The south neighbor of any cell in the first latitude band and the north neighbor of any cell in the last latitude band is `0` (sentinel value). Always check for `0` before using a neighbor ID as an array index.

4. **`edge_outward_normal` on degenerate edges.** For zero-length polar edges, the returned normal is `zero(SVector{3,Float64})` (the zero vector). The return type is always a `NamedTuple{(:base_point, :normal)}`, not a plain vector.

5. **No boundary.** Because the grid covers the full closed sphere, [`boundary_nodes`](@ref) and [`boundary_edges`](@ref) return empty vectors for any marker argument. There is no partial-sphere or regional mode.

6. **Internal copies.** The constructor copies the `lat_edges` and `lon_edges` arrays. Mutating the original arrays after construction has no effect on the grid.
