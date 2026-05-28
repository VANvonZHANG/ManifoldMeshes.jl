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
