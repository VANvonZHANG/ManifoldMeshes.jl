```@meta
CurrentModule = ManifoldMeshes
```

# ManifoldMeshes

Mesh infrastructure for scientific computing on manifolds. Currently provides latitude-longitude grids on the unit sphere (S²).

## Installation

```julia
using Pkg
Pkg.add("https://github.com/VANvonZHANG/ManifoldMeshes.jl")
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

# Cell geometry
cell_volume(grid, 1)
cell_centroid(grid, 1)

# Topological connectivity
cell_nodes(grid, 1)
cell_cells(grid, 1)
```

## Design Philosophy

**Let geometry belong to manifolds, topology belong to meshes.**

- All continuous mathematics (geodesics, distances, projections) is delegated to [Manifolds.jl](https://github.com/JuliaManifolds/Manifolds.jl)
- This library defines discrete connectivity (cells, nodes, edges) and computes geometric measures
- Cell boundaries are geodesic arcs (great circles), not coordinate lines

## API Reference

```@index
```

```@autodocs
Modules = [ManifoldMeshes]
```
