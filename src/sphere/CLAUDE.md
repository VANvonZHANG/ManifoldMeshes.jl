# sphere/ — Sphere Grid Implementations

Four S² grid types sharing spherical geometry utilities from `utils.jl`.

## Shared Utilities (`utils.jl`)

- `spherical_triangle_area(R, A, B, C)` — l'Huilier's formula, with `clamp(acos, -1, 1)` and `max(tan_half, 0.0)` guards
- All grid types compute `cell_volume` by splitting quadrilaterals into two spherical triangles

## LatLonGrid

Structured latitude-longitude grid, global coverage `[-90°, 90°] × [0°, 360°]`. Trait: `IsGrid`.

**Data layout:**
- `nodes::Matrix{SVector{3,Float64}}` — `(nlat+1) × (nlon+1)`, nodes[1, :] = south pole
- `cell_volumes::Matrix{Float64}` — `nlat × nlon`, pre-computed at construction
- `cell_centroids::Matrix{SVector{3,Float64}}` — `nlat × nlon`, pre-computed at construction
- Linear index: `cell_id = (ilat-1) * nlon + ilon`
- Node index: `node_id = (ilat-1) * (nlon+1) + ilon`
- Edge numbering: horizontal 1..(nlat+1)×nlon, vertical (nlat+1)×nlon+1..total

**Polar degeneration handling:**
- Each spherical quadrilateral is split along A-C diagonal into 2 spherical triangles
- If A-C degenerates, fall back to B-D diagonal
- If both diagonals degenerate (180° polar cell), use lune formula: `R² × (sin(lat₂) - sin(lat₁)) × |Δlon|`
- Zero-length edges return `zero(SVector{3,Float64})`

**Periodic boundaries:**
- `ilon_next = (ilon % nlon) + 1`, east-west wrapping
- North-south boundaries return 0 (sentinel)

## CubedSphereGrid

Cubed-sphere grid via gnomonic/equiangular projection mapping 6 faces onto S². Traits: `IsSemiGrid`, `MultiPatch`.

**Data layout:**
- `nodes::Vector{SVector{3,Float64}}` — node coordinates, shared nodes merged across face boundaries
- `n::Int` — cells per face per direction (total cells = 6n²)
- `_dual::RefValue` — lazy dual mesh cache

**Constructor:** `CubedSphereGrid(; n=4, projection=:gnomonic, R=1.0)`
- `projection`: `:gnomonic` (default) or `:equiangular`
- Faces numbered 1-6, in-face cells `(i, j)` mapped to global ID via `_cubed_sphere_cell_id(n, face, i, j)`

**Key functions:**
- `cell_face(g, cell_id)` → face ID (1-6)
- `cell_local_2d(g, cell_id)` → in-face `(i, j)` indices
- `_gnomonic_point(face, s, t)` / `_equiangular_point(face, s, t)` — face-local coords to S² mapping

**Cross-face topology:** `cell_cells` returns neighboring face's cell IDs at face boundaries, not 0

## ReducedGaussianGrid

Gaussian latitude band grid with longitude count decreasing toward poles for approximate equal-area cells. Trait: `IsSemiGrid`.

**Data layout:**
- `nodes::Vector{SVector{3,Float64}}` — all nodes (including poles)
- `_cell_nodes::Vector{NTuple{4,Int}}` — 4 node IDs per cell
- `_cell_edges::Vector{NTuple{4,Int}}` — 4 edge IDs per cell
- `_edge_nodes::Vector{NTuple{2,Int}}` — 2 endpoint IDs per edge
- `_cell_cells::Vector{NTuple{4,Int}}` — 4 neighbors (padded with 0 if fewer)
- `_dual::RefValue` — lazy dual mesh cache

**Constructor:** `ReducedGaussianGrid(; nlat=42, R=1.0)`
- `_gaussian_latitudes(n)` — computes Legendre roots as latitude boundaries
- `_octahedral_lon_counts(nlat)` — octahedral-style longitude reduction pattern

**Pole handling:** Polar cells degenerate into triangles (represented as quadrilaterals via duplicated nodes)

## HEALPixGrid

HEALPix Hierarchical Equal Area iso-Latitude Pixelization grid. Trait: `IsSemiGrid`.

**Data layout:**
- `nodes::Vector{SVector{3,Float64}}` — corner coordinates
- `nside::Int` — resolution parameter (total cells = 12 × nside²)
- `_cell_nodes::Vector{NTuple{4,Int}}` — with corner deduplication
- `_cell_edges::Vector{NTuple{4,Int}}`, `_cell_cells::Vector{NTuple{4,Int}}`
- `_edge_nodes::Vector{NTuple{2,Int}}`
- `_ordering::Symbol` — `:ring` or `:nested`
- `_dual::RefValue`

**Constructor:** `HEALPixGrid(; nside=4, ordering=:ring, R=1.0)`

**Morton/Z-order encoding:**
- `_morton_encode(ix, iy)` / `_morton_decode(m)` — interleaved bit coordinates
- `_nested_cell_center(nside, base, ix, iy)` — cell center coords in nested mode
- `_nested_to_ring(nside, nested_idx)` — nested→ring conversion (single source of truth)
- `_ring_to_nested_perm` field — cached ring→nested permutation (O(1) locate for `:nested`)

**Corner deduplication:** Shared corners between adjacent cells are merged to single node IDs at construction time

## Manifolds.jl Dependencies

- `Manifolds.mean(manifold, vertices)` — cell centroid (Fréchet mean)
- `Manifolds.mid_point(manifold, n1, n2)` — edge midpoint
- `Manifolds.project(manifold, base_point, vector)` — edge outward normal projected to tangent space
