# ManifoldMeshes.jl

Mesh infrastructure for scientific computing on manifolds. Currently provides latitude-longitude grids on the unit sphere (S²).

## Commands

```bash
# Test (from package root)
julia --project=. -e 'using Pkg; Pkg.test()'

# Run single test file
julia --project=. -e 'include("test/test_latlon_geometry.jl")'

# Aqua code quality checks
julia --project=. -e 'using Pkg; Pkg.test()'  # Aqua runs first in runtests.jl
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
├── traits.jl            # TopologyStyle (IsGrid/IsMesh), AbstractLocation (NodeLoc/CellLoc/EdgeLoc)
├── interface.jl         # AbstractManifoldMesh + 15 function stubs
└── sphere/
    └── latlon.jl        # LatLonGrid implementation
```

**Type hierarchy:**
- `AbstractManifoldMesh{M}` — parameterized by Manifolds.jl manifold type
- `LatLonGrid{M} <: AbstractManifoldMesh{M}` — only concrete implementation

**Interface functions** (all take `AbstractManifoldMesh`):
- Properties: `manifold`, `num_cells`, `num_nodes`, `num_edges`
- Geometry: `node_coordinates`, `cell_volume`, `cell_centroid`, `edge_length`, `edge_midpoint`, `edge_outward_normal`
- Topology: `cell_nodes`, `cell_cells`, `node_cells`, `cell_edges`
- Boundary: `boundary_nodes`, `boundary_edges`

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
- Test files: `test_traits.jl`, `test_latlon_{construction,geometry,connectivity,normals,edge_cases}.jl`, `test_performance.jl`

## Dependencies

- `Manifolds.jl` / `ManifoldsBase.jl` for all manifold geometry
- `StaticArrays.jl` for `SVector` and `NTuple` (zero GC pressure)
- Never add a new dependency without strong justification

## Gotchas

- Full-sphere grids have **no boundary** — `boundary_nodes` and `boundary_edges` return empty vectors
- Node at lon=0 and lon=360 are the **same physical point** with different linear IDs
- `edge_outward_normal` returns `NamedTuple{:base_point, :normal}` (tangent space semantics), not a plain vector

## Commit Convention

Use conventional commits: `feat:`, `fix:`, `chore:`, `docs:`, `ci:`, `test:`, `refactor:`
