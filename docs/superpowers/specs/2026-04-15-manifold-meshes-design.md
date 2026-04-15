# ManifoldMeshes.jl Design Specification

**Date**: 2026-04-15
**Status**: Approved
**Phase**: 1 (Sphere LatLon only, single-threaded, no I/O)

---

## 1. Design Philosophy

**Let geometry belong to manifolds, topology belong to meshes.**

- **Geometry**: All continuous mathematics (geodesics, distances, exponential/logarithmic maps, projections) is delegated to Manifolds.jl. This library never re-implements a geodesic.
- **Topology**: This library defines how discrete entities (cells, nodes, edges) connect on a manifold, and computes discrete geometric measures (cell volume, edge normal) from those connections.
- **No physics**: Physical quantities (temperature, velocity, emissions) are stored in a separate library (ManifoldFields.jl). A mesh struct never holds field data.
- **Pure manifold school**: Cell boundaries are geodesic arcs (great circles on the sphere), not coordinate lines. A "lat-lon cell" is a spherical quadrilateral whose edges are great-circle arcs connecting its four corner vertices.

## 2. Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Scope | Sphere S² only | Grid definition on arbitrary manifolds is an open problem; S² has unique properties (constant curvature, closed-form geodesics, exact area formulas) |
| Cell geometry | Geodesic polygons (great-circle edges) | Consistent with Manifolds.jl's distance/exp/log; volume via l'Huilier's formula |
| Return types for small fixed-size data | `SVector{N,Float64}`, `NTuple{K,Int}` | Stack allocation, zero GC pressure in tight PDE loops |
| Coordinate convention | Geographic latitude (equator=0, poles=±90), longitude (0..360) | Standard in geosciences |
| Caching | Pre-compute volumes, centroids, node positions at construction time | Space-for-time tradeoff; critical for FVM assembly performance |
| Traits vs inheritance | Traits (TopologyStyle, AbstractLocation) | Julian idiom; enables algorithm dispatch without deep type hierarchies |
| I/O | Not in Phase 1 | Keep the kernel minimal; downstream libraries or users handle file formats |

## 3. Dependencies

```
ManifoldsBase.jl  ◄── Manifolds.jl
       │                  │
       └──── ManifoldMeshes.jl  (this library)
```

| Package | Purpose |
|---------|---------|
| Manifolds.jl | Sphere type, distance, exp, log, mean, mid_point, project, is_point |
| ManifoldsBase.jl | AbstractManifold interface |
| StaticArrays.jl | SVector for stack-allocated coordinates |
| LinearAlgebra | cross product, dot product, normalize |

## 4. Type System

### 4.1 TopologyStyle Trait

Distinguishes structured grids (implicit topology) from unstructured meshes (explicit topology).

```julia
abstract type TopologyStyle end
struct IsGrid <: TopologyStyle end
struct IsMesh <: TopologyStyle end

# Type-level dispatch
TopologyStyle(::Type{<:AbstractManifoldMesh}) = IsMesh()

# Instance-level convenience (follows Base.IndexStyle convention)
TopologyStyle(m::AbstractManifoldMesh) = TopologyStyle(typeof(m))
```

### 4.2 AbstractLocation Trait

Tags for where physical data "lives" on a mesh. Defined here (not in ManifoldFields.jl) because boundary markers, connectivity queries, and node/edge/cell counts are mesh-level concepts.

```julia
abstract type AbstractLocation end
struct NodeLoc <: AbstractLocation end
struct CellLoc <: AbstractLocation end
struct EdgeLoc <: AbstractLocation end
```

### 4.3 AbstractManifoldMesh

```julia
abstract type AbstractManifoldMesh{M<:AbstractManifold} end
```

Parameterized by the underlying Manifolds.jl manifold type (e.g., `Sphere{2}`).

## 5. Interface Contract

Every mesh type must implement all functions below. No default implementations are provided — each mesh type defines its own optimized version.

### 5.1 Global Information

```julia
manifold(m::AbstractManifoldMesh)          -> AbstractManifold
num_cells(m::AbstractManifoldMesh)         -> Int
num_nodes(m::AbstractManifoldMesh)         -> Int
num_edges(m::AbstractManifoldMesh)         -> Int
```

### 5.2 Topological Connectivity

```julia
# Fixed-size returns (stack allocated)
cell_nodes(m, cell_id)  -> NTuple{N,Int}    # N vertices per cell
cell_cells(m, cell_id)  -> NTuple{K,Int}    # K neighbor cells
cell_edges(m, cell_id)  -> NTuple{N,Int}    # N edges per cell

# Variable-length returns (heap allocated — a node is shared by varying numbers of cells)
node_cells(m, node_id)  -> Vector{Int}
```

All indices are 1-based.

### 5.3 Geometry and Measures

All geometric computations delegate to Manifolds.jl primitives where possible.

```julia
# Coordinates: fixed SVector (stack allocated)
node_coordinates(m, node_id)  -> SVector{D,Float64}   # D = embedding dimension
cell_centroid(m, cell_id)    -> SVector{D,Float64}

# Scalar measures
cell_volume(m, cell_id)      -> Float64               # in steradians (or manifold-native units)
edge_length(m, edge_id)      -> Float64               # geodesic arc length

# Edge midpoint (independent query)
edge_midpoint(m, edge_id)    -> SVector{D,Float64}

# Outward unit normal: returns (base_point, tangent_vector)
# The tangent vector lives in T_{base_point}M and points away from the specified cell.
# Required for FVM flux computation.
edge_outward_normal(m, edge_id, cell_id) -> NamedTuple{(:base_point, :normal), Tuple{SVector{D,Float64}, SVector{D,Float64}}}
```

### 5.4 Boundary Markers

Markers can be `Symbol` (e.g., `:land`, `:ocean`) or `Int` (e.g., Gmsh physical group IDs). Boundaryless manifolds (e.g., full-sphere LatLon) return empty arrays.

```julia
boundary_nodes(m, marker)  -> Vector{Int}
boundary_edges(m, marker)  -> Vector{Int}
```

## 6. LatLonGrid Implementation

### 6.1 Struct Definition

```julia
struct LatLonGrid <: AbstractManifoldMesh{Sphere{2}}
    manifold::Sphere{2}
    lat_edges::Vector{Float64}    # latitude boundaries in degrees, ascending, [-90, ..., 90]
    lon_edges::Vector{Float64}    # longitude boundaries in degrees, ascending, [0, ..., 360]
    R::Float64                    # sphere radius (default 1.0)

    # Cached at construction time
    nlat::Int
    nlon::Int
    nodes::Matrix{SVector{3,Float64}}       # [ilat, ilon], size (nlat+1) × (nlon+1)
    cell_volumes::Matrix{Float64}            # [ilat, ilon], size nlat × nlon
    cell_centroids::Matrix{SVector{3,Float64}}  # [ilat, ilon], size nlat × nlon
end
```

### 6.2 Constructor

```julia
function LatLonGrid(; lat_edges::Vector{Float64}, lon_edges::Vector{Float64}, R::Float64=1.0)
    # Validate: lat_edges ascending from -90 to 90, lon_edges ascending from 0 to 360
    # Build nodes matrix via geographic→embedding conversion
    # Pre-compute cell_volumes via l'Huilier's formula (2 triangles per cell)
    # Pre-compute cell_centroids via Manifolds.mean on 4 vertices
    # Return struct
end
```

### 6.3 Cell Indexing Convention

Cell `(ilat, ilon)` has linear index `(ilat - 1) * nlon + ilon`.
- `ilat ∈ 1:nlat` (south to north)
- `ilon ∈ 1:nlon` (west to east)
- Cell 1 corresponds to the southwesternmost cell (nearest to lat=-90, lon=0)

Node `(ilat, ilon)` has linear index `(ilat - 1) * (nlon + 1) + ilon`.
- `ilat ∈ 1:(nlat+1)`, `ilon ∈ 1:(nlon+1)`

### 6.4 Coordinate Conversion

Geographic latitude → ℝ³ embedding (geocentric, Z = north pole):

```julia
function node_coordinates(g::LatLonGrid, node_id::Int)
    lat, lon = _node_latlon(g, node_id)
    θ, φ = deg2rad(lat), deg2rad(lon)
    return g.R * SVector(cos(θ)*cos(φ), cos(θ)*sin(φ), sin(θ))
end
```

### 6.5 Cell Volume (Pure Manifold School)

Cells are geodesic quadrilaterals. Volume computed by splitting into 2 spherical triangles and applying l'Huilier's formula:

```julia
function _spherical_triangle_area(R::Float64, A::SVector{3}, B::SVector{3}, C::SVector{3})
    a = acos(clamp(dot(B, C), -1, 1))
    b = acos(clamp(dot(A, C), -1, 1))
    c = acos(clamp(dot(A, B), -1, 1))
    s = (a + b + c) / 2
    tan_half = tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2)
    E = 4 * atan(sqrt(max(tan_half, 0.0)))  # spherical excess
    return R^2 * E
end
```

For each cell with vertices A, B, C, D (counterclockwise when viewed from outside):
```
area = triangle_area(A, B, C) + triangle_area(A, C, D)
```

### 6.6 Edge Length

```julia
function edge_length(g::LatLonGrid, edge_id::Int)
    n1, n2 = _edge_endpoints(g, edge_id)
    return Manifolds.distance(g.manifold, n1, n2)  # great-circle arc length
end
```

### 6.7 Edge Outward Normal

The outward normal to a great-circle arc at its midpoint:

1. **Polar collapse guard**: At the poles, multiple nodes map to the same ℝ³ point (e.g., all nodes at lat=90° become `(0,0,R)`). A degenerate edge with coincident endpoints has `cross(n1,n2) = 0`, causing NaN. Short-circuit: return zero normal for degenerate edges (flux through a zero-length edge is zero).
2. Compute midpoint via `Manifolds.mid_point(M, n1, n2)`
3. Great circle normal = `cross(n1, n2)` (perpendicular to the great circle plane)
4. Tangent along the arc = `normalize(cross(great_circle_normal, midpoint))`
5. Determine which side the cell is on: `sign(dot(great_circle_normal, cell_centroid))`
6. Outward normal = `cell_side * cross(tangent_along_arc, midpoint)`, projected to tangent space

```julia
function edge_outward_normal(g::LatLonGrid, edge_id::Int, cell_id::Int)
    n1, n2 = _edge_endpoints(g, edge_id)

    # Guard: polar collapse — degenerate edge with coincident endpoints
    if Manifolds.distance(g.manifold, n1, n2) < 1e-14
        midpoint = n1  # n1 ≈ n2
        return (base_point = SVector{3,Float64}(midpoint),
                normal = zero(SVector{3,Float64}))
    end

    midpoint = Manifolds.mid_point(g.manifold, n1, n2)
    gc_normal = cross(SVector(n1), SVector(n2))  # great circle plane normal

    tangent = normalize(cross(gc_normal, midpoint))
    cell_c = cell_centroid(g, cell_id)
    cell_side = sign(dot(gc_normal, SVector(cell_c)))

    outward = cell_side * cross(tangent, SVector(midpoint))
    outward = Manifolds.project(g.manifold, midpoint, outward)

    return (base_point = SVector{3,Float64}(midpoint),
            normal = SVector{3,Float64}(outward))
end
```

### 6.8 Connectivity (Implicit Topology for IsGrid)

All connectivity functions return fixed-size `NTuple` to ensure type stability. Missing neighbors (at poles) use `0` as sentinel; callers must check `if neighbor == 0 continue end`.

```julia
function cell_nodes(g::LatLonGrid, cell_id::Int)
    ilat, ilon = _cell_indices(g, cell_id)
    return (
        _node_linear_index(g, ilat,   ilon),     # SW
        _node_linear_index(g, ilat,   ilon + 1), # SE
        _node_linear_index(g, ilat+1, ilon + 1), # NE
        _node_linear_index(g, ilat+1, ilon),     # NW
    )
end

function cell_cells(g::LatLonGrid, cell_id::Int)
    ilat, ilon = _cell_indices(g, cell_id)
    south = ilat > 1     ? _cell_linear_index(g, ilat - 1, ilon) : 0
    north = ilat < nlat  ? _cell_linear_index(g, ilat + 1, ilon) : 0
    west  = _cell_linear_index(g, ilat, ilon == 1    ? nlon   : ilon - 1)  # periodic
    east  = _cell_linear_index(g, ilat, ilon == nlon ? 1      : ilon + 1)  # periodic
    return (south, north, west, east)  # always NTuple{4, Int}
end
```

Longitude wraps periodically. The grid has no boundaries (full sphere). Note: `cell_volume` handles degenerate polar cells correctly via l'Huilier's formula — when two vertices coincide (both at the pole), the degenerate triangle contributes zero area.

### 6.9 Boundary Markers

```julia
boundary_nodes(g::LatLonGrid, marker) = Int[]
boundary_edges(g::LatLonGrid, marker) = Int[]
```

A full-sphere LatLonGrid has no boundary. Polar singularities are handled implicitly through the geodesic quadrilateral geometry — no special marker needed.

### 6.10 TopologyStyle

```julia
TopologyStyle(::Type{LatLonGrid}) = IsGrid()
```

## 7. File Structure

```
ManifoldMeshes.jl/
├── Project.toml
├── src/
│   ├── ManifoldMeshes.jl       # entry point: using deps, include files, reexport
│   ├── traits.jl               # TopologyStyle, AbstractLocation
│   ├── interface.jl            # AbstractManifoldMesh + 15 interface function signatures
│   └── sphere/
│       └── latlon.jl           # LatLonGrid: struct, constructor, all interface implementations
└── test/
    ├── test_traits.jl
    ├── test_latlon_construction.jl
    ├── test_latlon_geometry.jl       # volume (lobe test), distance, centroid
    ├── test_latlon_connectivity.jl   # neighbor queries, periodicity
    └── test_latlon_normals.jl        # edge_outward_normal correctness
```

## 8. Tests

### 8.1 Lobe Test (Cell Volume Correctness)

Sum of all cell volumes must equal the total sphere surface area:
```
Σ cell_volume(g, i) == 4πR²
```

This is the most fundamental sanity check.

### 8.2 Symmetry Test

For a uniform LatLon grid (e.g., 90×180), all cells at the same latitude band must have equal volume.

### 8.3 Connectivity Periodicity

- Cell at `(ilat, 1)` must have a western neighbor at `(ilat, nlon)`.
- Cell at `(ilat, nlon)` must have an eastern neighbor at `(ilat, 1)`.

### 8.4 Normal Orthogonality

For each edge normal `(p, n)`:
- `dot(p, n) ≈ 0` (normal is in tangent space)
- `dot(n, n) ≈ 1` (unit normal)
- `dot(n, edge_tangent) ≈ 0` (normal perpendicular to edge)

### 8.5 Volume vs Geographic Formula Comparison

For a coarse grid near the equator, the geodesic quadrilateral area should closely match the geographic formula `R²·(sin lat₂ - sin lat₁)·Δlon`. Difference increases toward the poles.

### 8.6 Polar Collapse Robustness

- `edge_outward_normal` for a degenerate polar edge (zero-length) must return `normal = zero(SVector{3,Float64})` without NaN or error.
- `cell_volume` for polar cells (with two coincident pole vertices) must return the correct spherical triangle area (not zero).

## 9. Out of Scope (Phase 1)

- I/O (NetCDF, UGRID, Gmsh)
- Parallelism (threads, MPI)
- CubedSphereGrid, IcosahedralMesh, HEALPix
- Non-spherical manifolds
- Adaptive mesh refinement
- ManifoldFields.jl integration (separate library)
- ManifoldRegrid.jl integration (separate library)
- Vector field parallel transport
