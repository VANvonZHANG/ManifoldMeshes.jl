# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.5.0] - 2026-07-10

### Added

- **CSR topology storage** (`CSRMapping`): Compressed Sparse Row data structure for O(1) topology lookups, replacing O(N) linear scans across all 4 grid types
- **Edge topology**: `edge_cells`, `node_edges`, `edge_nodes` interface functions implemented for LatLonGrid, CubedSphereGrid, ReducedGaussianGrid, and HEALPixGrid
- **Batch query API**: Zero-copy batch functions (`all_cell_centroids`, `all_edge_lengths`, `all_cell_volumes`) returning `Tuple`-of-`SVector` for zero GC pressure
- **`node_cells` precomputation**: Inverse mapping (node→cells) precomputed at construction time for all grid types with `@inbounds` optimizations
- **Point location**: `locate_cell(g, lat, lon)` and 3D `locate_cell(g, p::SVector{3})` across LatLonGrid, CubedSphereGrid, ReducedGaussianGrid, HEALPixGrid
- **Interpolation weights**: `interpolation_weights(g, cell_id, lat, lon)` returning self-contained `(nodes, weights)` tuples for NodeLoc bilinear interpolation
- **Coordinate helpers**: internal `_latlon_to_cartesian` / `_cartesian_to_latlon` (shared with future `node_lonlat` UGRID accessor)

### Changed

- Migrated all 4 grid types to `CSRMapping`-based topology storage (LatLonGrid, CubedSphereGrid, ReducedGaussianGrid, HEALPixGrid)
- Extracted `_permute_uniform_csr` shared helper for reuse across uniform-cell grid types
- CSR accessors now return `SubArray` via `@view` for zero-copy semantics
- ReducedGaussianGrid: added `band_cell_offsets` and `node_lat_points` fields for O(1) point location (backwards-compatible; constructor-grown)

### Fixed

- Corrected LatLonGrid edge counts
- Fixed `_check_edge_id` placement in HEALPixGrid
- Optimized HEALPixGrid nested ordering permutation
- Eliminated `unique(cn)` allocation in ReducedGaussianGrid `node_cells` build
- Updated test type checks from `Vector{Int}` to `AbstractVector{Int}` for `@view` compatibility

## [0.4.0] - 2026-05-28

### Added

- **Generic visualization across all grid types**: `edge_segments()` and `cell_polygons()` now derive edges from `cell_nodes` boundary order, eliminating hardcoded quad-cell assumptions and working correctly for any grid topology
- **`cell_triangles()` function**: Decomposes each cell into K triangles fanning from the sphere-projected centroid; returns shared-vertex `(vertices, faces)` suitable for `GeometryBasics.Mesh`
- **Batched rendering**: Edge wireframes use a single NaN-separated `lines!()` call; filled meshes use a single `mesh!(verts, faces)` call with per-vertex coloring for improved performance
- **Extended visualization tests**: Both `test_visualization_data.jl` and `test_visualization_smoke.jl` now cover all 4 grid types; restored `slerp` and `node_points` unit tests
- **`examples/grids_visualization.jl`**: Standalone script producing a 4×2 comparison figure of all S² grid types

### Fixed

- Explicit imports and type stability in visualization module
- Hardened color array type for cross-backend compatibility

### Documentation

- **Visualization guide updates**: Added `cell_triangles()` API documentation with GLMakie/GeometryBasics usage example, batched rendering implementation note, and `examples/grids_visualization.jl` cross-reference
- **Tutorial and landing page**: Added Further Exploration subsection in Stage 5 and Visualization Examples link in Next Steps

## [0.3.0] - 2026-05-21

### Added

- **CubedSphereGrid**: cubed-sphere grid via gnomonic/equiangular projection with 6 faces, node merging across face boundaries, `cell_face` and `cell_local_2d` patch queries
- **ReducedGaussianGrid**: Gaussian latitude bands with pole-aware equal-area cell construction and reduced longitude count per band
- **HEALPixGrid**: HEALPix hierarchical equal-area pixelization with nested (Morton/Z-order) and ring ordering support
- **Trait system**: `IsSemiGrid` topology style, `CellTypeStyle` (`IsUniform{K}`/`IsMixed{MAX_K}`), `PatchStyle` (`NoPatch`/`MultiPatch`), `MixedCellTopology`
- **Dual mesh framework**: `AbstractDualMesh{M}`, lazy `dual()` computation with `RefValue` caching, `has_dual` query
- **Full topology for all grid types**: `cell_nodes`, `cell_cells`, `node_cells`, `cell_edges` implemented across all 4 grid types
- **Edge geometry for all grid types**: `edge_length`, `edge_midpoint`, `edge_outward_normal` for LatLonGrid, CubedSphereGrid, ReducedGaussianGrid, HEALPixGrid
- **Shared utilities**: `spherical_triangle_area` and `lune_area` extracted to `sphere/utils.jl`
- **Docstrings** for all grid constructors, patch queries, and trait types
- **Tests**: 387 tests covering all 4 grid types, traits, dual mesh, visualization data extraction, and performance baselines

### Changed

- Consolidated exhaustive test loops to spot-check pattern across all grid types (reduced from ~1700 to 387 tests while maintaining coverage)

## [0.2.0] - 2026-05-10

### Added

- 3D mesh visualization with `plot_mesh` and `plot_mesh_filled`, including semi-transparent sphere background rendering.
- `CairoMakie` as a runtime dependency for plotting support.
- Visualization smoke tests for CI coverage of plotting code.
- `VisualizationPlotting` module docstrings included in documentation manual.

### Changed

- Pinned CI and Documentation workflows to Julia 1.10 (LTS) for stability.
- TagBot workflow now includes daily cron schedule and CHANGELOG release notes extraction.

### Fixed

- TagBot permissions for automatic tag creation.
- Doc Preview Cleanup workflow permissions to prevent 403 errors.
- JuliaFormatter CI format check alignment.
- Bumped `CairoMakie` 0.12 → 0.15 and `GeometryBasics` 0.4 → 0.5 to fix Julia 1.12 AutoMerge precompilation failure.

## [0.1.1] - 2026-05-08

### Added

- Mesh data extraction layer for Makie/GeometryBasics visualization (`cell_polygons`, `edge_segments`, `node_markers`)
- CHANGELOG.md for tracking version history
- TagBot workflow integration with CHANGELOG release notes extraction

### Fixed

- `slerp` function vector comprehension bug in mesh visualization
- `cell_polygons` implementation for correct polygon vertex ordering

### Changed

- Applied JuliaFormatter with `sciml` style across the entire codebase
- Simplified CI to a single Ubuntu + Julia 1.11 runner for faster feedback
- Refined release workflow: removed standalone `release.yml`, relying on TagBot for automated GitHub releases

## [0.1.0] - 2025-04-19

### Added

- Initial release of ManifoldMeshes.jl
- `LatLonGrid` implementation for unit sphere (S²) with latitude-longitude discretization
- Full `AbstractManifoldMesh` interface: `num_cells`, `num_nodes`, `num_edges`, `node_coordinates`, `cell_volume`, `cell_centroid`, `edge_length`, `edge_midpoint`, `edge_outward_normal`
- Topology functions: `cell_nodes`, `cell_cells`, `node_cells`, `cell_edges`, `boundary_nodes`, `boundary_edges`
- Comprehensive test suite covering polar cells, periodic boundaries, and degenerate diagonals
- Documenter.jl documentation site with API reference
- GitHub Actions CI, CompatHelper, TagBot, and Dependabot automation

[unreleased]: https://github.com/VANvonZHANG/ManifoldMeshes.jl/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/VANvonZHANG/ManifoldMeshes.jl/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/VANvonZHANG/ManifoldMeshes.jl/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/VANvonZHANG/ManifoldMeshes.jl/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/VANvonZHANG/ManifoldMeshes.jl/releases/tag/v0.2.0
[0.1.1]: https://github.com/VANvonZHANG/ManifoldMeshes.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/VANvonZHANG/ManifoldMeshes.jl/releases/tag/v0.1.0
