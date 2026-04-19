# ManifoldMeshes.jl

Mesh infrastructure for scientific computing on manifolds. Currently provides latitude-longitude grids on the unit sphere (S²).

## Design Principles

- **Geometry belongs to manifolds, topology belongs to meshes.** All continuous mathematics (geodesics, distances, projections) is delegated to Manifolds.jl. This library never re-implements a geodesic.
- **No physics.** Physical quantities are stored in a separate library (ManifoldFields.jl). A mesh struct never holds field data.
- **Pure manifold school.** Cell boundaries are geodesic arcs (great circles), not coordinate lines.
- **Performance through caching.** Pre-compute volumes, centroids, and node positions at construction time.

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

## Dependencies

- `Manifolds.jl` / `ManifoldsBase.jl` for all manifold geometry
- `StaticArrays.jl` for `SVector` and `NTuple` (zero GC pressure)
- Never add a new dependency without strong justification

## Commit Convention

Use conventional commits: `feat:`, `fix:`, `chore:`, `docs:`, `ci:`, `test:`, `refactor:`
