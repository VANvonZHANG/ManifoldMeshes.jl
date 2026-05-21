```@meta
CurrentModule = ManifoldMeshes
```

# ManifoldMeshes.jl

Mesh infrastructure for scientific computing on Riemannian manifolds. Built on [Manifolds.jl](https://github.com/JuliaManifolds/Manifolds.jl) for all differential geometry — this library never re-implements a geodesic.

ManifoldMeshes provides discrete mesh structures (cells, nodes, edges) on the unit sphere **S²**, computing geometric measures from connections while delegating all continuous mathematics to the manifold.

## Grid Types

| Type | Description | Topology | Cells |
|------|-------------|----------|-------|
| [LatLonGrid](grids/latlon.md) | Structured latitude-longitude grid | `IsGrid` | Uniform quads |
| [CubedSphereGrid](grids/cubed_sphere.md) | Cubed-sphere via gnomonic projection, 6 faces | `IsSemiGrid`, `MultiPatch` | Uniform quads |
| [ReducedGaussianGrid](grids/reduced_gaussian.md) | Gaussian latitude bands, reduced longitude count | `IsSemiGrid` | Uniform quads |
| [HEALPixGrid](grids/healpix.md) | Hierarchical equal-area iso-latitude pixels | `IsSemiGrid` | Uniform quads |

## Quick Start

```julia
using ManifoldMeshes

# Create a latitude-longitude grid (5×10 cells)
grid = LatLonGrid(lat_edges=range(-90, 90; length=6), lon_edges=range(0, 360; length=11))

# Query cell geometry
cell_volume(grid, 1)    # spherical area of cell 1
cell_centroid(grid, 1)  # SVector{3,Float64} on the sphere

# Explore connectivity
cell_nodes(grid, 1)     # node IDs bounding cell 1
cell_cells(grid, 1)     # neighboring cell IDs

# Switch grid type — same interface
hp = HEALPixGrid(nside=4)
cell_volume(hp, 1)
cell_nodes(hp, 1)
```

## Next Steps

- **[Tutorial](tutorial.md)** — step-by-step walkthrough from grid creation to visualization
- **[Theory](theory.md)** — why geodesic cells, how grids are constructed, and the math behind the measures
- **[Grid Types](grids/latlon.md)** — detailed guide for each grid type
- **[API Reference](api.md)** — complete function and type documentation
