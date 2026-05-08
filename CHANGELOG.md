# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Visualization data extraction layer: `slerp`, `node_points`, `edge_segments`, `cell_polygons` for Makie ecosystem integration.
- `GeometryBasics` dependency (compat 0.4) for `Point3f` interop.

### Fixed

- `slerp` vector comprehension bug that returned 3n points instead of n.
- `cell_polygons` reverse arc concatenation logic to produce closed polygons.

### Changed

- Simplified CI to a single runner (Julia 1.11 + ubuntu-latest), removing macOS and multi-version matrix.
- Downgraded `GeometryBasics` compat from 0.5 to 0.4 for Makie ecosystem compatibility.

## [0.1.0] - 2026-05-08

### Added

- `LatLonGrid` — latitude-longitude grid on the unit sphere with full manifold mesh interface.
- Geometry queries: `node_coordinates`, `cell_volume`, `cell_centroid`, `edge_length`, `edge_midpoint`, `edge_outward_normal`.
- Topology queries: `cell_nodes`, `cell_cells`, `node_cells`, `cell_edges`.
- Boundary markers: `boundary_nodes`, `boundary_edges` (empty for full-sphere grids).
- Pre-computed cached volumes, centroids, and node positions at construction time.
- Comprehensive test suite covering construction, geometry, connectivity, normals, edge cases, and performance.
- Aqua code quality checks integrated into test suite.
- Documenter.jl documentation with GitHub Pages deployment.
- GitHub Actions CI, CompatHelper, TagBot, and Dependabot workflows.
- JuliaFormatter configuration (sciml style).

[unreleased]: https://github.com/VANvonZHANG/ManifoldMeshes.jl/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/VANvonZHANG/ManifoldMeshes.jl/releases/tag/v0.1.0
