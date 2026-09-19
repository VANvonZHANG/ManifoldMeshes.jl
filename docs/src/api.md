```@meta
CurrentModule = ManifoldMeshes
```

# API Reference

## Types

### Mesh Types

```@docs
AbstractManifoldMesh
LatLonGrid
CubedSphereGrid
ReducedGaussianGrid
HEALPixGrid
AbstractDualMesh
```

### Traits

```@docs
TopologyStyle
IsGrid
IsSemiGrid
IsMesh
CellTypeStyle
IsUniform
IsMixed
PatchStyle
NoPatch
MultiPatch
MixedCellTopology
```

### Locations

```@docs
AbstractLocation
NodeLoc
CellLoc
EdgeLoc
```

## Properties

```@docs
manifold
num_cells
num_nodes
num_edges
```

## Geometry

```@docs
node_coordinates
cell_volume
cell_centroid
```

## Spherical Polygon Geometry

Predicates, areas, and clipping for the geodesic polygons that bound cells.
Rings are unit vectors in `cell_nodes` order, counter-clockwise seen from
outside the sphere; `cell_ring` builds them from a mesh.

```@docs
side_of_geodesic
geodesic_arc_intersection
cell_ring
spherical_triangle_area
spherical_polygon_area
spherical_polygon_intersection
```

## Topology

```@docs
cell_nodes
cell_cells
node_cells
cell_edges
```

## Edge Geometry

```@docs
edge_length
edge_midpoint
edge_outward_normal
```

## Boundary

```@docs
boundary_nodes
boundary_edges
```

## Dual Mesh

```@docs
dual
has_dual
```

## Patch System

```@docs
cell_face
cell_local_2d
```

## Visualization

```@docs
slerp
node_points
edge_segments
cell_polygons
cell_triangles
plot_mesh
plot_mesh_filled
```
