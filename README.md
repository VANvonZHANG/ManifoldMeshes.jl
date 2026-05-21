# ManifoldMeshes

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://VANvonZHANG.github.io/ManifoldMeshes.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://VANvonZHANG.github.io/ManifoldMeshes.jl/dev/)
[![Build Status](https://github.com/VANvonZHANG/ManifoldMeshes.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/VANvonZHANG/ManifoldMeshes.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/VANvonZHANG/ManifoldMeshes.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/VANvonZHANG/ManifoldMeshes.jl)

Mesh infrastructure for scientific computing on Riemannian manifolds. Built on [Manifolds.jl](https://github.com/JuliaManifolds/Manifolds.jl) — this library never re-implements a geodesic.

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

grid = LatLonGrid(lat_edges=range(-90, 90; length=6), lon_edges=range(0, 360; length=11))
cell_volume(grid, 1)    # spherical area
cell_centroid(grid, 1)  # SVector{3,Float64}
cell_nodes(grid, 1)     # node IDs bounding cell 1
```

**[Full Documentation](https://VANvonZHANG.github.io/ManifoldMeshes.jl/dev/)** — tutorial, theory, per-grid guides, and complete API reference.

## License

MIT
