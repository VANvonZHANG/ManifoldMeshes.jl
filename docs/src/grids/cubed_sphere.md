# CubedSphereGrid

`CubedSphereGrid` maps the six faces of an inscribed cube onto the unit sphere $S^2$ via a gnomonic (or equiangular) central projection.
Each face is a structured $n \times n$ grid of spherical quadrilaterals, yielding $6n^2$ cells total with no polar singularity.
All node positions, cell volumes, and cell centroids are pre-computed at construction time, making every subsequent query an O(1) lookup.

```@raw html
<table style="margin: 1em 0;">
<tr><th style="text-align:left; padding-right:1.5em;">Trait</th><th style="text-align:left;">Value</th></tr>
<tr><td style="padding-right:1.5em;"><code>TopologyStyle</code></td><td><code>IsGrid()</code></td></tr>
<tr><td style="padding-right:1.5em;"><code>CellTypeStyle</code></td><td><code>IsUniform{4}()</code></td></tr>
<tr><td style="padding-right:1.5em;"><code>PatchStyle</code></td><td><code>MultiPatch{6}()</code></td></tr>
</table>
```

## Construction

```julia
CubedSphereGrid(; n, projection=:gnomonic, rotation=I, R=1.0)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `n` | `Int` | Number of cells per face edge ($\geq 1$). Total cells = $6n^2$. |
| `projection` | `Symbol` | `:gnomonic` (default) or `:equiangular`. Controls how the cube faces are projected onto the sphere. |
| `rotation` | `SMatrix{3,3,Float64,9}` | Rotation matrix applied to all nodes (default: identity). Useful for orienting the grid relative to a geographic or geophysical frame. |
| `R` | `Float64` | Sphere radius (default `1.0`). Must be positive. |

```@example cubed
using ManifoldMeshes

# n=3: 3 cells per face edge, 6 faces => 54 cells total
grid = CubedSphereGrid(n = 3)

println("Cells: ", num_cells(grid))    # 54
println("Nodes: ", num_nodes(grid))    # 96
println("Edges: ", num_edges(grid))    # 144
```

Using the equiangular projection, which distributes area more evenly across each face:

```@example cubed
grid_eq = CubedSphereGrid(n = 4, projection = :equiangular)

println("Cells: ", num_cells(grid_eq))  # 96
```

Rotating the entire grid by 90 degrees around the Z axis:

```@example cubed
using StaticArrays

rot = SMatrix{3,3}(0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
grid_rot = CubedSphereGrid(n = 2, rotation = rot)

println("Cells: ", num_cells(grid_rot)) # 24
```

## Cell Geometry

### Shape

Each cell is a spherical quadrilateral whose four vertices are obtained by projecting face-local coordinates $(s, t) \in [-1, 1]^2$ onto $S^2$.
For the gnomonic projection, a point on face $f$ at local coordinates $(s, t)$ is mapped to

$$\mathbf{p}(s, t) = R \cdot \frac{\mathbf{v}_f(s, t)}{\|\mathbf{v}_f(s, t)\|}$$

where $\mathbf{v}_f$ is the un-normalized cube-face vector (e.g., for face 1: $\mathbf{v}_1 = (s, t, 1)$).
The equiangular projection replaces $(s, t)$ with $(\tan^{-1} s, \tan^{-1} t)$ before mapping, producing a more uniform area distribution.

Cell nodes are returned in the order **SW, SE, NE, NW** (southwest, southeast, northeast, northwest) within the face-local coordinate system.

### Volume

Cell volumes (solid angles on the unit sphere, scaled by $R^2$) are pre-computed at construction using l'Huilier's formula.
Each quadrilateral is split into two spherical triangles along the SW-NE diagonal, and their areas are summed.
The gnomonic projection introduces slight area variation: cells near the center of a face are slightly smaller than cells near the edges.

```@example cubed
g = CubedSphereGrid(n = 4)

# Volume of cell 1 (near the center of face 1)
v = cell_volume(g, 1)
println("Volume of cell 1: ", round(v; digits=6))

# Compare a center cell vs. an edge cell on face 1 (n=4)
# Cell (i=3, j=3) on face 1 => cell_id = (3-1)*4 + 3 = 11
# Cell (i=4, j=1) on face 1 => cell_id = (1-1)*4 + 4 = 4
v_center = cell_volume(g, 11)
v_edge = cell_volume(g, 4)
println("Center cell: ", round(v_center; digits=6))
println("Edge cell:   ", round(v_edge; digits=6))
```

### Centroid

Cell centroids are the Fréchet mean (intrinsic mean) of the four vertices on $S^2$, computed via `Manifolds.mean`.
They are returned as `SVector{3, Float64}` in Cartesian coordinates and always lie on the sphere surface.

```@example cubed
c = cell_centroid(g, 1)
println("Centroid of cell 1: ", round.(c; digits=4))
```

### Area Conservation

The sum of all cell volumes equals $4\pi R^2$ to machine precision, regardless of projection type or resolution.

```@example cubed
total = sum(cell_volume(g, i) for i in 1:num_cells(g))
println("Sum of volumes: ", round(total; digits=6))
println("4*pi*R^2:       ", round(4pi * g.R^2; digits=6))
```

## Topology and Connectivity

Cells are numbered sequentially by face: cells $1$ through $n^2$ belong to face 1, cells $n^2 + 1$ through $2n^2$ to face 2, and so on.
Within each face, cell indices follow row-major order:

$$\text{cell\_id} = (\text{face} - 1) \cdot n^2 + (j - 1) \cdot n + i$$

### Cell-to-Cell Neighbors

[`cell_cells`](@ref) returns a 4-tuple `(south, north, west, east)` of neighbor cell IDs.
For interior cells (not on a face boundary), all four neighbors are valid cell IDs within the same face.
At face boundaries, the missing cross-face neighbor is reported as `0` (sentinel value).

```@example cubed
g4 = CubedSphereGrid(n = 4)

# Interior cell on face 3 (i=2, j=2)
cell_id = (3 - 1) * 16 + (2 - 1) * 4 + 2  # cell 26
south, north, west, east = cell_cells(g4, cell_id)
println("Interior neighbors of cell $cell_id: S=$south, N=$north, W=$west, E=$east")

# Cell at the east edge of face 1 (i=4, j=2)
edge_cell = (2 - 1) * 4 + 4  # cell 8
neighbors = cell_cells(g4, edge_cell)
println("Edge cell $edge_cell neighbors: ", neighbors)
println("Has boundary (0): ", any(n -> n == 0, neighbors))
```

### Cell-to-Node Mapping

[`cell_nodes`](@ref) returns the four corner node IDs as `(SW, SE, NE, NW)`.
Each face owns its own set of $(n+1)^2$ nodes, so nodes on shared edges and corners have different IDs even though they coincide geometrically.

```@example cubed
nodes = cell_nodes(g4, 1)
println("Nodes of cell 1: ", nodes)
```

### Cell-to-Edge Mapping

[`cell_edges`](@ref) returns four edge IDs as `(south, north, west, east)`.
Edges are numbered per face: horizontal edges first ($(n+1) \times n$ per face), then vertical edges.

```@example cubed
edges = cell_edges(g4, 1)
println("Edges of cell 1: ", edges)
```

### Boundary

Because CubedSphereGrid covers the full sphere, [`boundary_nodes`](@ref) and [`boundary_edges`](@ref) always return empty vectors for any marker.

```@example cubed
println("Boundary nodes: ", boundary_nodes(g4, :default))
println("Boundary edges: ", boundary_edges(g4, :default))
```

## Patch System

CubedSphereGrid is the only grid type with the `MultiPatch{6}` trait, meaning it is composed of 6 independent structured patches (the cube faces).
Two patch-specific query functions are available:

### `cell_face`

[`cell_face`](@ref) returns the face index (1--6) containing a given cell.

```@example cubed
g = CubedSphereGrid(n = 3)

# Face assignment is sequential: cells 1..9 => face 1, 10..18 => face 2, etc.
println("Cell 1 => face ", cell_face(g, 1))     # face 1
println("Cell 10 => face ", cell_face(g, 10))   # face 2
println("Cell 54 => face ", cell_face(g, 54))   # face 6

# Verify: each face has exactly n^2 cells
for f in 1:6
    n_cells = count(i -> cell_face(g, i) == f, 1:num_cells(g))
    println("Face $f: $n_cells cells")
end
```

### `cell_local_2d`

[`cell_local_2d`](@ref) returns the within-face 2D indices `(i, j)` of a cell, where both `i` and `j` range from 1 to `n`.

```@example cubed
g = CubedSphereGrid(n = 3)

# First cell on face 1: local (1,1)
println("Cell 1 local: ", cell_local_2d(g, 1))

# Last cell on face 1 (n=3): local (3,3)
println("Cell 9 local: ", cell_local_2d(g, 9))

# First cell on face 2: local (1,1)
println("Cell 10 local: ", cell_local_2d(g, 10))
```

These two functions enable efficient patch-local algorithms: you can extract an $n \times n$ subgrid from a single face, run structured numerical methods on it, then stitch results back together.

## When to Use

### Strengths

- **No polar singularity.** Unlike [LatLonGrid](latlon.md), no coordinate lines converge anywhere on the sphere. All cells remain well-conditioned quadrilaterals, making the cubed sphere suitable for global atmospheric and ocean models that require uniform CFL conditions.
- **Quasi-uniform cell area.** The gnomonic projection distributes area much more evenly than latitude-longitude grids. The equiangular projection further reduces area variation.
- **Structured per face.** Within each face, cells have regular $(i, j)$ indexing, enabling efficient stencil-based numerical methods (finite differences, finite volumes).
- **Patch system.** [`cell_face`](@ref) and [`cell_local_2d`](@ref) provide natural support for distributed-memory parallelism where each MPI rank owns one or more faces.
- **Conservation.** The sum of all cell volumes equals $4\pi R^2$ to machine precision.

### Weaknesses

- **Cross-face connectivity.** Cells at face boundaries report `0` for their cross-face neighbors in [`cell_cells`](@ref), rather than the adjacent face's cell ID. Algorithms that traverse neighbor chains must handle the `0` sentinel explicitly.
- **Duplicate boundary nodes.** Nodes on shared edges and corners are not deduplicated: coincident points have different IDs. This means topological queries across faces require additional bookkeeping.
- **Slight area variation.** The gnomonic projection does not produce perfectly equal-area cells. If exact equal area is needed, consider [HEALPixGrid](healpix.md) instead.
- **Non-trivial coordinate system.** The face-local $(s, t)$ coordinates and the projection mapping are more complex than simple $(\theta, \phi)$ indexing.

### Typical Use Cases

- Global atmospheric circulation models (e.g., CAM-FV, GEOScubed)
- Ocean general circulation models requiring quasi-uniform resolution
- Parallel simulations where the 6-face patch structure maps naturally to domain decomposition
- Any application needing uniform resolution without polar singularities

## Gotchas

1. **Cross-face neighbors are `0`.** [`cell_cells`](@ref) returns `0` for neighbors that would lie on an adjacent face, rather than computing the actual cross-face cell ID. Always check for `0` before using a neighbor ID as an array index.

2. **Nodes are not deduplicated across face boundaries.** Each face stores its own copy of edge and corner nodes. Two nodes from adjacent faces that coincide geometrically will have different IDs. This affects topology: the number of nodes is $6(n+1)^2$, not the smaller count you would get with shared boundary nodes.

3. **Face ordering is fixed.** Face normals are: 1: +Z, 2: -Z, 3: +Y, 4: -Y, 5: +X, 6: -X. This determines the orientation and placement of each face on the sphere.

4. **`node_cells` is O(cells).** The implementation scans all cells to find those adjacent to a given node. For performance-sensitive code with frequent node-to-cell queries, pre-compute a lookup table.

5. **No boundary.** Because the grid covers the full closed sphere, [`boundary_nodes`](@ref) and [`boundary_edges`](@ref) return empty vectors for any marker argument.

6. **`edge_outward_normal` return type.** The returned value is a `NamedTuple{(:base_point, :normal)}`, not a plain vector. The `base_point` is the edge midpoint on the sphere and `normal` is the outward-pointing tangent-space vector.

7. **Gnomonic vs. equiangular projection.** The default `:gnomonic` projection has greater area distortion near face edges. The `:equiangular` projection reduces this distortion by applying $\arctan$ to the local coordinates. Both projections satisfy area conservation globally ($\sum V_i = 4\pi R^2$), but individual cell areas differ between the two.
