# UnstructuredMesh

```@docs
UnstructuredMesh
```

`UnstructuredMesh` stores an arbitrary spherical polygon mesh as *data*: node
coordinates plus a padded face-node connectivity table. The four parametric
grid types generate all geometry from constructor parameters;
`UnstructuredMesh` instead derives it from the table, which is what makes
foreign meshes (MPAS, ICON) first-class.

## Construction

```julia
using ManifoldMeshes

# From a UGRID-style 0-based table; pad inactive slots with fill_value
mesh = UnstructuredMesh(node_lon, node_lat, face_nodes; R = 1.0,
                        start_index = 0, fill_value = -1)

# From Cartesian points on the unit sphere
mesh = UnstructuredMesh(Sphere(2), points, face_nodes)
```

Corners must be listed in cyclic boundary order (starting corner irrelevant);
cells may be triangles, quads, pentagons, hexagons, ... (per-cell arity,
`IsMixed{MAX_K}`). v1 supports `Sphere` manifolds only, and cells must be
convex.

## Point location

`locate_cell` builds a k-d tree over cell centroids lazily (first call) and
refines with an exact spherical point-in-polygon test — O(log n) per query vs
the analytic O(1) of parametric grids. Points on shared boundaries resolve to
the smallest containing cell id.

## Interpolation

`interpolation_weights` is arity-aware: quad cells use bilinear weights from
local coordinates solved on the cell's own corners (gnomonic projection +
Newton — the same machinery as `HEALPixGrid`), other arities use Wachspress
coordinates (which reduce to barycentric coordinates on triangles).
