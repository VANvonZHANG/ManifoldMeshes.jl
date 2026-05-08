# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-05-08

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
