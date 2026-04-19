# ManifoldMeshes

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://VANvonZHANG.github.io/ManifoldMeshes.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://VANvonZHANG.github.io/ManifoldMeshes.jl/dev/)
[![Build Status](https://github.com/VANvonZHANG/ManifoldMeshes.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/VANvonZHANG/ManifoldMeshes.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/VANvonZHANG/ManifoldMeshes.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/VANvonZHANG/ManifoldMeshes.jl)

Mesh infrastructure for scientific computing on manifolds. Currently provides latitude-longitude grids on the unit sphere (S²).

## Design Philosophy

**Let geometry belong to manifolds, topology belong to meshes.**

- All continuous mathematics (geodesics, distances, projections) is delegated to [Manifolds.jl](https://github.com/JuliaManifolds/Manifolds.jl)
- This library defines discrete connectivity (cells, nodes, edges) and computes geometric measures (cell volume, edge normals) from those connections
- Cell boundaries are geodesic arcs (great circles), not coordinate lines

## Installation

```julia
] add https://github.com/VANvonZHANG/ManifoldMeshes.jl
```

## Quick Start

```julia
using ManifoldMeshes

# Create a 5×10 latitude-longitude grid on the unit sphere
grid = LatLonGrid(lat_edges=range(-90, 90; length=6), lon_edges=range(0, 360; length=11))

# Global information
num_cells(grid)   # 50
num_nodes(grid)   # 66
num_edges(grid)   # 120

# Query a cell (1-indexed, ilat × nlon + ilon)
cell_id = 1
cell_volume(grid, cell_id)       # spherical area of the cell
cell_centroid(grid, cell_id)     # centroid as SVector{3,Float64}
cell_nodes(grid, cell_id)        # (1, 2, 13, 12)
cell_cells(grid, cell_id)        # neighboring cell IDs (periodic in longitude)

# Node coordinates
node_coordinates(grid, 1)        # South pole SVector{3,Float64}

# Edge geometry
edge_id = 1
edge_length(grid, edge_id)       # great-circle arc length
edge_midpoint(grid, edge_id)     # midpoint on sphere
```

## Current Features

- `LatLonGrid` on `Sphere(2)` with configurable resolution
- Pre-computed cell volumes (spherical triangle decomposition via l'Huilier's formula)
- Pre-computed node positions and cell centroids (3D Cartesian)
- Topological connectivity: cell-cell, cell-node, node-cell, cell-edges
- Edge operations: length, midpoint, outward normal (tangent-space projected)
- Periodic longitude boundary (seamless wrapping)
- Polar collapse handling (degenerate edges/normals at poles)

## License

MIT
