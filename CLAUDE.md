# ManifoldMeshes.jl

Mesh infrastructure for scientific computing on manifolds. Provides structured and semi-structured grids on the unit sphere (S²): `LatLonGrid`, `CubedSphereGrid`, `ReducedGaussianGrid`, and `HEALPixGrid`.

## Commands

```bash
# Test (from package root)
julia --project=. -e 'using Pkg; Pkg.test()'

# Run single test file
julia --project=. -e 'include("test/test_latlon_geometry.jl")'

# Aqua code quality checks
julia --project=. -e 'using Pkg; Pkg.test()'  # Aqua runs first in runtests.jl

# Visualization (requires adding a Makie backend to your environment)
using CairoMakie  # or GLMakie for interactive 3D
fig = plot_mesh(grid; view=Symbol("3d"))
save("mesh.png", fig)
```

## Design Principles

- **Geometry belongs to manifolds, topology belongs to meshes.** All continuous mathematics (geodesics, distances, projections) is delegated to Manifolds.jl. This library never re-implements a geodesic.
- **No physics.** Physical quantities are stored in a separate library (ManifoldFields.jl). A mesh struct never holds field data.
- **Pure manifold school.** Cell boundaries are geodesic arcs (great circles), not coordinate lines.
- **Performance through caching.** Pre-compute volumes, centroids, and node positions at construction time.

## Architecture

```
src/
├── ManifoldMeshes.jl    # Module entry, exports, includes
├── traits.jl            # TopologyStyle, CellTypeStyle, PatchStyle, AbstractLocation, MixedCellTopology
├── interface.jl         # AbstractManifoldMesh + full function interface
├── dual.jl              # AbstractDualMesh + lazy dual mesh framework
├── sphere/
│   ├── utils.jl         # Shared: _spherical_triangle_area, _lune_area
│   ├── latlon.jl        # LatLonGrid (structured lat-lon, IsGrid)
│   ├── cubed_sphere.jl  # CubedSphereGrid (gnomonic projection, 6-face, IsSemiGrid)
│   ├── reduced_gaussian.jl  # ReducedGaussianGrid (Gaussian lat bands, IsSemiGrid)
│   └── healpix.jl       # HEALPixGrid (Nside hierarchical, IsSemiGrid)
└── visualization/
    ├── mesh_data.jl     # Data extraction: node_points, edge_segments, cell_polygons
    └── plotting.jl      # plot_mesh, plot_mesh_filled (requires Makie at runtime)
```

**Type hierarchy:**
- `AbstractManifoldMesh{M}` — parameterized by Manifolds.jl manifold type
  - `LatLonGrid{M} <: AbstractManifoldMesh{M}` — structured lat-lon (`IsGrid`)
  - `CubedSphereGrid{M} <: AbstractManifoldMesh{M}` — cubed-sphere (`IsSemiGrid`, `MultiPatch`)
  - `ReducedGaussianGrid{M} <: AbstractManifoldMesh{M}` — Gaussian lat bands (`IsSemiGrid`)
  - `HEALPixGrid{M} <: AbstractManifoldMesh{M}` — HEALPix hierarchical (`IsSemiGrid`)
  - `AbstractDualMesh{M} <: AbstractManifoldMesh{M}` — lazy dual mesh marker

**Trait system:**
- `TopologyStyle`: `IsGrid` | `IsSemiGrid` | `IsMesh`
- `CellTypeStyle`: `IsUniform{K}` | `IsMixed{MAX_K}`
- `PatchStyle`: `NoPatch` | `MultiPatch`
- `AbstractLocation`: `NodeLoc` | `CellLoc` | `EdgeLoc`

**Interface functions** (all take `AbstractManifoldMesh`):
- Properties: `manifold`, `num_cells`, `num_nodes`, `num_edges`
- Geometry: `node_coordinates`, `cell_volume`, `cell_centroid`
- Topology: `cell_nodes`, `cell_cells`, `node_cells`, `cell_edges`
- Edge: `edge_length`, `edge_midpoint`, `edge_outward_normal`
- Boundary: `boundary_nodes`, `boundary_edges`
- Dual: `dual`, `has_dual`
- Patch (cubed-sphere): `cell_face`, `cell_local_2d`
- Visualization: `slerp`, `node_points`, `edge_segments`, `cell_polygons`, `plot_mesh`, `plot_mesh_filled`

## Code Style

- Formatter: `.JuliaFormatter.toml` with `style = "sciml"` (4-space indent, aligned assignments)
- Functions: snake_case (`cell_volume`, `edge_outward_normal`)
- Types: PascalCase (`LatLonGrid`, `AbstractManifoldMesh`)
- Internal helpers: prefix with `_` (`_spherical_triangle_area`, `_cell_linear_index`)
- Return types for small fixed-size data: `SVector{N,Float64}`, `NTuple{K,Int}`

## Testing Requirements

- Every public function must have a test
- Always test edge cases: polar cells, periodic boundaries, degenerate diagonals
- Fundamental sanity checks: `Σ cell_volume = 4πR²`
- Spot-check pattern: avoid exhaustive loops over all cells; test representative indices
- Test files: `test_traits.jl`, `test_latlon_{construction,geometry,connectivity,normals,edge_cases}.jl`, `test_cubed_sphere.jl`, `test_reduced_gaussian.jl`, `test_healpix.jl`, `test_dual.jl`, `test_performance.jl`, `test_visualization_{data,smoke}.jl`

## Dependencies

- `Manifolds.jl` / `ManifoldsBase.jl` for all manifold geometry
- `StaticArrays.jl` for `SVector` and `NTuple` (zero GC pressure)
- `CairoMakie.jl` / `GeometryBasics.jl` for visualization (should become weak deps)
- Never add a new dependency without strong justification

## Gotchas

- Full-sphere grids have **no boundary** — `boundary_nodes` and `boundary_edges` return empty vectors
- Node at lon=0 and lon=360 are the **same physical point** with different linear IDs (LatLonGrid)
- `edge_outward_normal` returns `NamedTuple{:base_point, :normal}` (tangent space semantics), not a plain vector
- `plot_mesh` requires Makie to be loaded **before** calling — run `using CairoMakie` or `using GLMakie` first
- HEALPix uses Morton/Z-order curve for nested ordering — cell IDs follow hierarchical pixel numbering
- CubedSphereGrid nodes are merged across face boundaries — shared edge nodes have single global IDs
- ReducedGaussianGrid has `IsUniform{4}` cells but variable node counts per latitude band
- `dual()` computes lazily on first call and caches the result in `_dual::RefValue`

## Commit Convention

Use conventional commits: `feat:`, `fix:`, `chore:`, `docs:`, `ci:`, `test:`, `refactor:`
