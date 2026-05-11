# TODO — ManifoldMeshes.jl Development Roadmap

> Last updated: 2026-05-11  
> Current version: v0.2.0  
> Status: Core S² lat-lon grid is solid; expansion toward general manifolds and unstructured meshes is the next frontier.

---

## 1. Manifold Coverage (New Mesh Types)

The library currently has **one** concrete implementation: `LatLonGrid` on S². The abstract interface is designed to support arbitrary manifolds.

| Priority | Item | Description |
|----------|------|-------------|
| 🔴 High | **Cubed-sphere grid on S²** | Equal-area cells, avoids polar singularity of lat-lon. Essential for climate/weather codes. |
| 🔴 High | **Icosahedral / hexagonal grid on S²** | Hex-pentagon dual mesh. Standard in geodesic dome / atmospheric modeling. |
| 🟡 Medium | **Toroidal grid (T²)** | Periodic in both dimensions. Test-bed for non-spherical manifolds. |
| 🟡 Medium | **Stereographic / projective plane grids** | Non-orientable manifold examples. |
| 🟢 Low | **Higher-dimensional spheres (S³, S⁴)** | Requires extending the interface to handle 3D/4D cells (tetrahedra, simplices). |

---

## 2. Topology & Mesh Types

Currently only `IsGrid` (structured) is implemented. The `IsMesh` trait exists but has no concrete subtype.

| Priority | Item | Description |
|----------|------|-------------|
| 🔴 High | **Unstructured triangular mesh on S²** | General polygon mesh (`IsMesh`). Each cell can have arbitrary node count. |
| 🔴 High | **Dual mesh construction** | Given a primal mesh, compute the dual (nodes ↔ cells, edges ↔ edges). Critical for staggered-grid schemes. |
| 🟡 Medium | **Mesh refinement (quad-tree / tri-tree)** | Adaptive local refinement with hanging-node handling. |
| 🟡 Medium | **Mesh coarsening / agglomeration** | Inverse of refinement; useful for multigrid. |
| 🟢 Low | **Mesh I/O (VTK, NetCDF, JSON)** | Read/write mesh topology and geometry to standard formats. |

---

## 3. Interface Gaps

The current 15-function interface is minimal. Several common mesh operations are missing.

| Priority | Item | Description |
|----------|------|-------------|
| 🔴 High | **Batch query functions** | `all_cell_volumes(g) -> Vector{Float64}`, `all_node_coordinates(g) -> Vector{SVector{3,Float64}}`, etc. Current per-cell allocation pattern is slow in tight loops. |
| 🔴 High | **`edge_cells(g, edge_id)`** | Return the 1 or 2 cells adjacent to an edge. Essential for flux computations in FVM. |
| 🟡 Medium | **`node_edges(g, node_id)`** | Return all edges incident to a node. |
| 🟡 Medium | **Diagonal neighbor queries** | `cell_cells` currently only returns face-sharing neighbors. Add option for vertex-sharing (diagonal) neighbors. |
| 🟡 Medium | **Proper boundary marker system** | `boundary_nodes(g, marker)` and `boundary_edges(g, marker)` currently ignore `marker` and always return `Int[]` for closed manifolds. For open manifolds (e.g. hemispheres), implement named boundary segments. |
| 🟢 Low | **Mesh validation utilities** | Check for overlapping cells, negative volumes, dangling nodes, inconsistent orientation. |
| 🟢 Low | **Mesh quality metrics** | Aspect ratio, skewness, orthogonality, max/min edge ratio per cell. |

---

## 4. Performance & Memory

| Priority | Item | Description |
|----------|------|-------------|
| 🔴 High | **CSR/CSG sparse layout for topology** | `node_cells` currently allocates a `Vector{Int}` per call. Pre-compute CSR arrays for allocation-free adjacency walks. |
| 🟡 Medium | **Lazy / on-demand geometry caching** | Currently all volumes, centroids, and node positions are computed at construction. For very fine grids, allow lazy evaluation with LRU cache. |
| 🟡 Medium | **Static array return types for batch queries** | When grid dimensions are known at compile time (e.g. `LatLonGrid{6,11}`), return `SVector` or `NTuple` instead of `Vector`. |
| 🟢 Low | **GPU / CUDA support** | Store node coordinates and cell data on GPU; dispatch geometry queries to CUDA kernels. |
| 🟢 Low | **Distributed mesh partitioning** | Domain decomposition for MPI parallelization (e.g. METIS/ParMETIS integration). |

---

## 5. Visualization Enhancements

| Priority | Item | Description |
|----------|------|-------------|
| 🟡 Medium | **Interactive 3D plots (GLMakie backend)** | Current `plot_mesh` works with CairoMakie. Add GLMakie interactivity (rotation, zoom, cell picking). |
| 🟡 Medium | **Color mapping by scalar field** | `plot_mesh_filled` currently supports per-cell flat color. Add continuous colormap based on vertex-interpolated or cell-centered scalar data. |
| 🟡 Medium | **Vector field visualization** | Overlay arrows / streamlines on the mesh surface for tangent vector fields. |
| 🟢 Low | **Animation support** | Time-evolving mesh plots with `Makie.Record`. |
| 🟢 Low | **Export to PNG/SVG/PDF** | Document-quality figure export helpers with proper DPI and LaTeX fonts. |

---

## 6. Dependencies & Packaging

| Priority | Item | Description |
|----------|------|-------------|
| 🔴 High | **Convert CairoMakie to a weak dependency (Package Extension)** | Currently `CairoMakie` is in `[deps]`, forcing all users to install it even if they only need mesh topology. Use Julia 1.9+ extensions so visualization is opt-in. |
| 🟡 Medium | **Compat bounds audit** | Ensure all `[compat]` entries are tight enough to prevent future breakages, especially for Manifolds.jl which is still evolving. |
| 🟢 Low | **Reduce runtime dependencies** | Evaluate whether `GeometryBasics` is needed outside visualization; if not, move it to the extension as well. |

---

## 7. Documentation

| Priority | Item | Description |
|----------|------|-------------|
| 🔴 High | **API reference with `@autodocs`** | Currently docs/src/ only has `index.md`. Add API pages for `AbstractManifoldMesh`, `LatLonGrid`, `VisualizationPlotting`, and all interface functions. |
| 🔴 High | **Tutorial: Building a custom mesh type** | Step-by-step guide implementing `MyMesh <: AbstractManifoldMesh` for a new manifold. |
| 🟡 Medium | **Comparison with other Julia mesh libraries** | Document differences vs. `Meshes.jl`, `Gridap.jl`, `Ferrite.jl`. |
| 🟡 Medium | **Performance best-practices guide** | When to use batch queries, how to avoid allocations in loops, memory layout tips. |
| 🟢 Low | **Visualization cookbook** | Gallery of common plot types with copy-paste code. |

---

## 8. Testing & Quality

| Priority | Item | Description |
|----------|------|-------------|
| 🟡 Medium | **Benchmark suite with BenchmarkTools** | Systematic benchmarking of `cell_volume`, `edge_outward_normal`, `node_cells`, etc. Track regressions in CI. |
| 🟡 Medium | **Property-based testing (Hypothesis.jl)** | Generate random valid `lat_edges` / `lon_edges` and assert invariants (total volume = 4πR², all normals point outward, etc.). |
| 🟢 Low | **Fuzz testing for invalid inputs** | Ensure graceful error messages for malformed constructor arguments. |

---

## 9. Ecosystem Integration

| Priority | Item | Description |
|----------|------|-------------|
| 🟡 Medium | **Interoperability with `GeometryBasics.Mesh`** | Provide `to_geometrybasics(g) -> GeometryBasics.Mesh` for downstream packages that consume mesh data. |
| 🟡 Medium | **Interoperability with `Meshes.jl`** | Conversion functions to/from `Meshes.jl` representations where applicable. |
| 🟢 Low | **Interoperability with `Trixi.jl` / `Oceananigans.jl`** | Adapter modules for popular PDE solvers. |

---

## Legend

- 🔴 **High**: Blocks core functionality or user adoption; tackle next.
- 🟡 **Medium**: Important for production use; schedule for next minor release.
- 🟢 **Low**: Nice-to-have; community contributions welcome.

---

## How to Contribute

1. Pick an item from the table above.
2. Open a draft PR referencing this TODO file.
3. Follow the existing code style (`.JuliaFormatter.toml`, `style = "sciml"`).
4. Every new public function needs a test in `test/`.
