# ManifoldMeshes

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://VANvonZHANG.github.io/ManifoldMeshes.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://VANvonZHANG.github.io/ManifoldMeshes.jl/dev/)
[![Build Status](https://github.com/VANvonZHANG/ManifoldMeshes.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/VANvonZHANG/ManifoldMeshes.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/VANvonZHANG/ManifoldMeshes.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/VANvonZHANG/ManifoldMeshes.jl)

Mesh infrastructure for scientific computing on Riemannian manifolds. Built on [Manifolds.jl](https://github.com/JuliaManifolds/Manifolds.jl) for all differential geometry — this library never re-implements a geodesic.

## Design Philosophy

**Let geometry belong to manifolds, topology belong to meshes.**

- All continuous mathematics (geodesics, distances, projections) is delegated to Manifolds.jl
- This library defines discrete connectivity (cells, nodes, edges) and computes geometric measures from those connections
- Cell boundaries are geodesic arcs (great circles), not coordinate lines

## Grid Types

| Type | Description | Cells | Traits |
|------|-------------|-------|--------|
| `LatLonGrid` | Structured latitude-longitude grid | Uniform quads | `IsGrid` |
| `CubedSphereGrid` | Cubed-sphere via gnomonic projection | Uniform quads, 6 faces | `IsSemiGrid`, `MultiPatch` |
| `ReducedGaussianGrid` | Gaussian latitude bands with reduced longitude count | Uniform quads | `IsSemiGrid` |
| `HEALPixGrid` | Hierarchical Equal Area iso-Latitude Pixelization | Uniform quads | `IsSemiGrid` |

## Installation

```julia
] add https://github.com/VANvonZHANG/ManifoldMeshes.jl
```

## Quick Start

```julia
using ManifoldMeshes

# LatLonGrid: 5×10 cells on the unit sphere
grid = LatLonGrid(lat_edges=range(-90, 90; length=6), lon_edges=range(0, 360; length=11))
num_cells(grid)   # 50

# CubedSphereGrid: 4×4 cells per face (96 total)
cs = CubedSphereGrid(n=4)

# ReducedGaussianGrid: T42 Gaussian latitudes
rg = ReducedGaussianGrid(nlat=42)

# HEALPixGrid: Nside=4 (192 cells)
hp = HEALPixGrid(nside=4)

# Query a cell
cell_id = 1
cell_volume(grid, cell_id)       # spherical area
cell_centroid(grid, cell_id)     # SVector{3,Float64}
cell_nodes(grid, cell_id)        # node IDs bounding the cell
cell_cells(grid, cell_id)        # neighboring cell IDs
```

## Interface

All grid types implement the `AbstractManifoldMesh` interface:

**Properties:** `manifold`, `num_cells`, `num_nodes`, `num_edges`

**Geometry:** `node_coordinates`, `cell_volume`, `cell_centroid`, `edge_length`, `edge_midpoint`, `edge_outward_normal`

**Topology:** `cell_nodes`, `cell_cells`, `node_cells`, `cell_edges`

**Boundary:** `boundary_nodes`, `boundary_edges` (empty for full-sphere grids)

**Dual mesh:** `dual(g)`, `has_dual(g)` — lazy construction

## Visualization

```julia
using CairoMakie  # or GLMakie for interactive 3D
fig = plot_mesh(grid; view=Symbol("3d"))
save("mesh.png", fig)
```

## Testing

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## License

MIT
