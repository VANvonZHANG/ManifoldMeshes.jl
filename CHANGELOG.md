# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[unreleased]: https://github.com/VANvonZHANG/ManifoldMeshes.jl/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/VANvonZHANG/ManifoldMeshes.jl/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/VANvonZHANG/ManifoldMeshes.jl/releases/tag/v0.2.0
[0.1.1]: https://github.com/VANvonZHANG/ManifoldMeshes.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/VANvonZHANG/ManifoldMeshes.jl/releases/tag/v0.1.0
