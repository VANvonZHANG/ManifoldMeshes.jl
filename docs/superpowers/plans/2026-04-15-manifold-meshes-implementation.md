# ManifoldMeshes.jl Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement Phase 1 of ManifoldMeshes.jl — a Julia library for sphere S² grid infrastructure with pure-manifold-school geodesic geometry, traits-based dispatch, and pre-computed caches for FVM performance.

**Architecture:** Traits (TopologyStyle, AbstractLocation) define dispatch categories. AbstractManifoldMesh{M} is the abstract supertype parameterized by Manifolds.jl manifold type. LatLonGrid is the sole concrete type in Phase 1: a structured grid on Sphere(2) with great-circle cell boundaries, pre-computed node positions, cell volumes (l'Huilier's formula), and cell centroids (Manifolds.mean).

**Tech Stack:** Julia 1.10, Manifolds.jl (Sphere, distance, mid_point, project, mean), StaticArrays.jl (SVector, NTuple), LinearAlgebra (cross, dot, normalize)

---

## File Structure

```
ManifoldMeshes.jl/
├── Project.toml                    # already exists
├── src/
│   ├── ManifoldMeshes.jl           # entry point: using deps, include files
│   ├── traits.jl                   # TopologyStyle, AbstractLocation, IsGrid, IsMesh
│   ├── interface.jl                # AbstractManifoldMesh{M}, docstrings for interface
│   └── sphere/
│       └── latlon.jl               # LatLonGrid: struct, constructor, all 15 interface functions
└── test/
    ├── runtests.jl                 # includes all test files
    ├── test_traits.jl
    ├── test_latlon_construction.jl # constructor, counts, node_coordinates
    ├── test_latlon_geometry.jl     # cell_volume (lobe test), cell_centroid
    ├── test_latlon_connectivity.jl # cell_nodes, cell_cells, node_cells, cell_edges
    └── test_latlon_normals.jl      # edge_length, edge_midpoint, edge_outward_normal, boundary
```

---

## Edge Numbering Scheme (LatLonGrid)

Total edges = `(nlat+1)*nlon + nlat*nlon`

| Type | Range | Grid indices | Connects |
|------|-------|-------------|----------|
| Horizontal | `1 : (nlat+1)*nlon` | `ilat ∈ 1:(nlat+1), ilon ∈ 1:nlon` | `node(ilat, ilon)` → `node(ilat, ilon%nlon+1)` |
| Vertical | `(nlat+1)*nlon+1 : end` | `ilat ∈ 1:nlat, ilon ∈ 1:nlon` | `node(ilat, ilon)` → `node(ilat+1, ilon)` |

Cell `(ilat, ilon)` edges (south, north, west, east):
- South: `(ilat-1)*nlon + ilon`
- North: `ilat*nlon + ilon`
- West: `(nlat+1)*nlon + (ilat-1)*nlon + ilon`
- East: `(nlat+1)*nlon + (ilat-1)*nlon + (ilon==nlon ? 1 : ilon+1)`

---

## Task 1: Project Skeleton + Entry Point

**Files:**
- Create: `src/ManifoldMeshes.jl`
- Create: `test/runtests.jl`

- [ ] **Step 1: Create the entry point**

```julia
# src/ManifoldMeshes.jl
module ManifoldMeshes

using Manifolds
using ManifoldsBase
using StaticArrays
using LinearAlgebra

include("traits.jl")
include("interface.jl")
include("sphere/latlon.jl")

end # module
```

- [ ] **Step 2: Create the test runner**

```julia
# test/runtests.jl
using ManifoldMeshes
using Test

@testset "ManifoldMeshes.jl" begin
    include("test_traits.jl")
    include("test_latlon_construction.jl")
    include("test_latlon_geometry.jl")
    include("test_latlon_connectivity.jl")
    include("test_latlon_normals.jl")
end
```

- [ ] **Step 3: Initialize git and install dependencies**

```bash
cd /home/zhangfan/Project/20260328_HEMCOManifold/ManifoldMeshes.jl
git init
git add Project.toml src/ManifoldMeshes.jl test/runtests.jl
git commit -m "init: project skeleton with entry point and test runner"
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Expected: `julia` resolves all dependencies. `Pkg.test()` will fail because no source files are included yet (traits.jl missing) — that's expected at this stage.

---

## Task 2: Traits (TopologyStyle, AbstractLocation)

**Files:**
- Create: `src/traits.jl`
- Create: `test/test_traits.jl`

- [ ] **Step 1: Write the failing test**

```julia
# test/test_traits.jl
@testset "TopologyStyle trait" begin
    @test TopologyStyle(IsGrid) === IsGrid()
    @test TopologyStyle(IsMesh) === IsMesh()

    # Instance-level forwarding should be defined (will work once LatLonGrid exists)
    # Tested here with direct type-level dispatch
    @test isa(TopologyStyle(LatLonGrid), TopologyStyle)  # type-level on abstract fallback
end

@testset "AbstractLocation types" begin
    @test NodeLoc <: AbstractLocation
    @test CellLoc <: AbstractLocation
    @test EdgeLoc <: AbstractLocation
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /home/zhangfan/Project/20260328_HEMCOManifold/ManifoldMeshes.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `TopologyStyle`, `IsGrid`, `IsMesh`, `AbstractLocation`, `NodeLoc`, `CellLoc`, `EdgeLoc`, `LatLonGrid` are all undefined.

- [ ] **Step 3: Implement traits**

```julia
# src/traits.jl
"""
    TopologyStyle

Trait distinguishing structured grids (implicit topology) from unstructured meshes (explicit topology).
"""
abstract type TopologyStyle end

"""Structured grid with implicit topology (e.g., LatLonGrid)."""
struct IsGrid <: TopologyStyle end

"""Unstructured mesh with explicit topology (not yet implemented)."""
struct IsMesh <: TopologyStyle end

# Type-level dispatch (follows Base.IndexStyle convention)
TopologyStyle(::Type{T}) where T = IsMesh()

# Instance-level forwarding
TopologyStyle(m::T) where T = TopologyStyle(T)

"""
    AbstractLocation

Trait tagging where physical data "lives" on a mesh.
"""
abstract type AbstractLocation end

"""Data at mesh nodes (vertices)."""
struct NodeLoc <: AbstractLocation end

"""Data at mesh cell centers."""
struct CellLoc <: AbstractLocation end

"""Data at mesh edges."""
struct EdgeLoc <: AbstractLocation end
```

- [ ] **Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for the `TopologyStyle trait` and `AbstractLocation types` testsets. The `LatLonGrid` line in the test will error since LatLonGrid doesn't exist yet — **fix the test**: remove or comment out the `TopologyStyle(LatLonGrid)` line (it will be re-tested in Task 4).

Updated test (remove the LatLonGrid line):

```julia
# test/test_traits.jl
@testset "TopologyStyle trait" begin
    @test TopologyStyle(IsGrid) === IsGrid()
    @test TopologyStyle(IsMesh) === IsMesh()
    @test isa(TopologyStyle(IsGrid), TopologyStyle)
end

@testset "AbstractLocation types" begin
    @test NodeLoc <: AbstractLocation
    @test CellLoc <: AbstractLocation
    @test EdgeLoc <: AbstractLocation
end
```

Re-run:
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/traits.jl test/test_traits.jl
git commit -m "feat: add TopologyStyle and AbstractLocation traits"
```

---

## Task 3: AbstractManifoldMesh + Interface

**Files:**
- Create: `src/interface.jl`

- [ ] **Step 1: Create the interface module**

No tests needed — this is purely abstract definitions with docstrings.

```julia
# src/interface.jl
"""
    AbstractManifoldMesh{M<:AbstractManifold}

Abstract supertype for all manifold meshes. Parameterized by the underlying
Manifolds.jl manifold type.

Every concrete subtype must implement the full interface defined below.
"""
abstract type AbstractManifoldMesh{M<:AbstractManifold} end

# ── Global Information ──────────────────────────────────────────────
# manifold(m) -> AbstractManifold
# num_cells(m) -> Int
# num_nodes(m) -> Int
# num_edges(m) -> Int

# ── Topological Connectivity ────────────────────────────────────────
# cell_nodes(m, cell_id)  -> NTuple{N,Int}
# cell_cells(m, cell_id)  -> NTuple{K,Int}
# cell_edges(m, cell_id)  -> NTuple{N,Int}
# node_cells(m, node_id)  -> Vector{Int}

# ── Geometry and Measures ───────────────────────────────────────────
# node_coordinates(m, node_id)  -> SVector{D,Float64}
# cell_centroid(m, cell_id)    -> SVector{D,Float64}
# cell_volume(m, cell_id)      -> Float64
# edge_length(m, edge_id)      -> Float64
# edge_midpoint(m, edge_id)    -> SVector{D,Float64}
# edge_outward_normal(m, edge_id, cell_id) -> NamedTuple{(:base_point, :normal)}

# ── Boundary Markers ───────────────────────────────────────────────
# boundary_nodes(m, marker)  -> Vector{Int}
# boundary_edges(m, marker)  -> Vector{Int}
```

- [ ] **Step 2: Commit**

```bash
git add src/interface.jl
git commit -m "feat: add AbstractManifoldMesh interface with docstrings"
```

---

## Task 4: LatLonGrid Struct + Constructor + Index Helpers + node_coordinates + Global Info

**Files:**
- Create: `src/sphere/latlon.jl`
- Modify: `test/test_traits.jl` (add LatLonGrid TopologyStyle test)
- Create: `test/test_latlon_construction.jl`

- [ ] **Step 1: Write the failing test**

```julia
# test/test_latlon_construction.jl
@testset "LatLonGrid construction" begin
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 90.0, 180.0, 270.0, 360.0])

    @test g.nlat == 2
    @test g.nlon == 4
    @test g.R == 1.0

    # Global info
    @test manifold(g) isa Sphere{2}
    @test num_cells(g) == 2 * 4
    @test num_nodes(g) == 3 * 5
    @test num_edges(g) == 3 * 4 + 2 * 4  # (nlat+1)*nlon + nlat*nlon

    # TopologyStyle
    @test TopologyStyle(g) === IsGrid()
    @test TopologyStyle(LatLonGrid) === IsGrid()
end

@testset "node_coordinates basics" begin
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 180.0, 360.0])

    # South pole: all nodes at lat=-90 map to (0, 0, -R)
    for ilon in 1:3
        n = node_coordinates(g, (1-1)*3 + ilon)
        @test n ≈ SVector(0.0, 0.0, -1.0) atol=1e-12
    end

    # North pole: all nodes at lat=90 map to (0, 0, R)
    for ilon in 1:3
        n = node_coordinates(g, (3-1)*3 + ilon)
        @test n ≈ SVector(0.0, 0.0, 1.0) atol=1e-12
    end

    # Equator at lon=0
    n = node_coordinates(g, (2-1)*3 + 1)  # node(2, 1): lat=0, lon=0
    @test n ≈ SVector(1.0, 0.0, 0.0) atol=1e-12

    # Equator at lon=90
    n = node_coordinates(g, (2-1)*3 + 2)  # node(2, 2): lat=0, lon=90
    @test n ≈ SVector(0.0, 1.0, 0.0) atol=1e-12

    # Equator at lon=180
    n = node_coordinates(g, (2-1)*3 + 3)  # node(2, 3): lat=0, lon=180
    @test n ≈ SVector(-1.0, 0.0, 0.0) atol=1e-12
end

@testset "constructor validation" begin
    # Missing -90 start
    @test_throws ArgumentError LatLonGrid(lat_edges=[-45.0, 0.0, 90.0], lon_edges=[0.0, 360.0])

    # Missing 90 end
    @test_throws ArgumentError LatLonGrid(lat_edges=[-90.0, 0.0, 45.0], lon_edges=[0.0, 360.0])

    # Not ascending
    @test_throws ArgumentError LatLonGrid(lat_edges=[-90.0, 90.0, 0.0], lon_edges=[0.0, 360.0])

    # Missing 0 start for lon
    @test_throws ArgumentError LatLonGrid(lat_edges=[-90.0, 90.0], lon_edges=[90.0, 360.0])

    # Missing 360 end for lon
    @test_throws ArgumentError LatLonGrid(lat_edges=[-90.0, 90.0], lon_edges=[0.0, 180.0])

    # Not ascending lon
    @test_throws ArgumentError LatLonGrid(lat_edges=[-90.0, 90.0], lon_edges=[0.0, 360.0, 180.0])

    # Need at least 2 edges each
    @test_throws ArgumentError LatLonGrid(lat_edges=[-90.0], lon_edges=[0.0, 360.0])
    @test_throws ArgumentError LatLonGrid(lat_edges=[-90.0, 90.0], lon_edges=[0.0])
end
```

Also add the TopologyStyle test for LatLonGrid to `test/test_traits.jl`:

```julia
# Add to the end of test/test_traits.jl, inside the TopologyStyle testset:
@testset "TopologyStyle for LatLonGrid" begin
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 360.0])
    @test TopologyStyle(g) === IsGrid()
    @test TopologyStyle(LatLonGrid) === IsGrid()
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `LatLonGrid`, `manifold`, `num_cells`, `num_nodes`, `num_edges`, `node_coordinates` undefined.

- [ ] **Step 3: Implement LatLonGrid**

```julia
# src/sphere/latlon.jl
using Manifolds: Sphere

"""
    LatLonGrid <: AbstractManifoldMesh{Sphere{2}}

Structured latitude-longitude grid on the sphere S² with geodesic (great-circle) cell boundaries.

# Fields
- `manifold::Sphere{2}` — the underlying manifold
- `lat_edges::Vector{Float64}` — latitude boundaries in degrees, ascending [-90, ..., 90]
- `lon_edges::Vector{Float64}` — longitude boundaries in degrees, ascending [0, ..., 360]
- `R::Float64` — sphere radius
- `nlat::Int` — number of latitude bands (cells in lat direction)
- `nlon::Int` — number of longitude sectors (cells in lon direction)
- `nodes::Matrix{SVector{3,Float64}}` — cached node positions, size (nlat+1) × (nlon+1)
- `cell_volumes::Matrix{Float64}` — cached cell volumes (steradians), size nlat × nlon
- `cell_centroids::Matrix{SVector{3,Float64}}` — cached cell centroids, size nlat × nlon
"""
struct LatLonGrid <: AbstractManifoldMesh{Sphere{2}}
    manifold::Sphere{2}
    lat_edges::Vector{Float64}
    lon_edges::Vector{Float64}
    R::Float64
    nlat::Int
    nlon::Int
    nodes::Matrix{SVector{3,Float64}}
    cell_volumes::Matrix{Float64}
    cell_centroids::Matrix{SVector{3,Float64}}
end

# ── Constructor ─────────────────────────────────────────────────────

function LatLonGrid(; lat_edges::Vector{Float64}, lon_edges::Vector{Float64}, R::Float64=1.0)
    # Validate lat_edges
    length(lat_edges) < 2 && throw(ArgumentError("lat_edges must have at least 2 elements"))
    length(lon_edges) < 2 && throw(ArgumentError("lon_edges must have at least 2 elements"))
    lat_edges[1] != -90.0 && throw(ArgumentError("lat_edges must start at -90"))
    lat_edges[end] != 90.0 && throw(ArgumentError("lat_edges must end at 90"))
    !issorted(lat_edges) && throw(ArgumentError("lat_edges must be ascending"))
    lon_edges[1] != 0.0 && throw(ArgumentError("lon_edges must start at 0"))
    lon_edges[end] != 360.0 && throw(ArgumentError("lon_edges must end at 360"))
    !issorted(lon_edges) && throw(ArgumentError("lon_edges must be ascending"))

    nlat = length(lat_edges) - 1
    nlon = length(lon_edges) - 1
    M = Sphere(2)

    # Build nodes matrix
    nodes = Matrix{SVector{3,Float64}}(undef, nlat + 1, nlon + 1)
    for ilat in 1:(nlat + 1)
        lat = lat_edges[ilat]
        θ = deg2rad(lat)
        sinθ, cosθ = sin(θ), cos(θ)
        for ilon in 1:(nlon + 1)
            lon = lon_edges[ilon]
            φ = deg2rad(lon)
            nodes[ilat, ilon] = R * SVector(cosθ * cos(φ), cosθ * sin(φ), sinθ)
        end
    end

    # Pre-compute cell volumes
    cell_volumes = Matrix{Float64}(undef, nlat, nlon)
    for ilat in 1:nlat
        A = nodes[ilat, 1]       # SW
        D = nodes[ilat + 1, 1]   # NW
        for ilon in 1:nlon
            ilon_next = ilon == nlon ? 1 : ilon + 1
            B = nodes[ilat, ilon_next]         # SE
            C = nodes[ilat + 1, ilon_next]     # NE
            cell_volumes[ilat, ilon] =
                _spherical_triangle_area(R, A, B, C) +
                _spherical_triangle_area(R, A, C, D)
        end
    end

    # Pre-compute cell centroids via Manifolds.mean
    cell_centroids = Matrix{SVector{3,Float64}}(undef, nlat, nlon)
    for ilat in 1:nlat
        for ilon in 1:nlon
            ilon_next = ilon == nlon ? 1 : ilon + 1
            verts = [nodes[ilat, ilon], nodes[ilat, ilon_next],
                     nodes[ilat+1, ilon_next], nodes[ilat+1, ilon]]
            c = Manifolds.mean(M, verts)
            cell_centroids[ilat, ilon] = SVector{3,Float64}(c)
        end
    end

    return LatLonGrid(M, lat_edges, lon_edges, R, nlat, nlon, nodes, cell_volumes, cell_centroids)
end

# ── Internal: Spherical Triangle Area (l'Huilier's formula) ────────

function _spherical_triangle_area(R::Float64, A::SVector{3,Float64},
                                  B::SVector{3,Float64}, C::SVector{3,Float64})
    a = acos(clamp(dot(B, C) / (R * R), -1, 1))
    b = acos(clamp(dot(A, C) / (R * R), -1, 1))
    c = acos(clamp(dot(A, B) / (R * R), -1, 1))
    s = (a + b + c) / 2
    tan_half = tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2)
    E = 4 * atan(sqrt(max(tan_half, 0.0)))  # spherical excess
    return R^2 * E
end

# ── Internal: Index Helpers ────────────────────────────────────────

function _cell_indices(g::LatLonGrid, cell_id::Int)
    ilat = div(cell_id - 1, g.nlon) + 1
    ilon = rem(cell_id - 1, g.nlon) + 1
    return (ilat, ilon)
end

function _cell_linear_index(g::LatLonGrid, ilat::Int, ilon::Int)
    return (ilat - 1) * g.nlon + ilon
end

function _node_indices(g::LatLonGrid, node_id::Int)
    ilat = div(node_id - 1, g.nlon + 1) + 1
    ilon = rem(node_id - 1, g.nlon + 1) + 1
    return (ilat, ilon)
end

function _node_linear_index(g::LatLonGrid, ilat::Int, ilon::Int)
    return (ilat - 1) * (g.nlon + 1) + ilon
end

# ── TopologyStyle Override ─────────────────────────────────────────

TopologyStyle(::Type{LatLonGrid}) = IsGrid()

# ── Global Information ─────────────────────────────────────────────

manifold(g::LatLonGrid) = g.manifold
num_cells(g::LatLonGrid) = g.nlat * g.nlon
num_nodes(g::LatLonGrid) = (g.nlat + 1) * (g.nlon + 1)
function num_edges(g::LatLonGrid)
    return (g.nlat + 1) * g.nlon + g.nlat * g.nlon
end

# ── Geometry: node_coordinates ─────────────────────────────────────

function node_coordinates(g::LatLonGrid, node_id::Int)
    ilat, ilon = _node_indices(g, node_id)
    return g.nodes[ilat, ilon]
end
```

- [ ] **Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add src/sphere/latlon.jl test/test_latlon_construction.jl test/test_traits.jl
git commit -m "feat: add LatLonGrid struct, constructor, node_coordinates, global info"
```

---

## Task 5: Cell Volume Correctness

**Files:**
- Create: `test/test_latlon_geometry.jl`

- [ ] **Step 1: Write the failing test**

```julia
# test/test_latlon_geometry.jl
@testset "spherical triangle area" begin
    # Octant of unit sphere: vertices at (1,0,0), (0,1,0), (0,0,1)
    A = SVector(1.0, 0.0, 0.0)
    B = SVector(0.0, 1.0, 0.0)
    C = SVector(0.0, 0.0, 1.0)
    area = ManifoldMeshes._spherical_triangle_area(1.0, A, B, C)
    @test area ≈ π / 2 atol=1e-12  # octant = π/2 steradians
end

@testset "lobe test: total area = 4πR²" begin
    for R in [1.0, 6371.0]
        g = LatLonGrid(lat_edges=collect(-90.0:10.0:90.0),
                       lon_edges=collect(0.0:15.0:360.0), R=R)
        total = sum(i -> cell_volume(g, i), 1:num_cells(g))
        @test total ≈ 4 * π * R^2 rtol=1e-10
    end
end

@testset "symmetry: same latitude band = same volume" begin
    g = LatLonGrid(lat_edges=collect(-90.0:10.0:90.0),
                   lon_edges=collect(0.0:15.0:360.0))
    for ilat in 1:g.nlat
        ilon1 = 1
        ilon2 = g.nlon
        id1 = _cell_linear_index(g, ilat, ilon1)
        id2 = _cell_linear_index(g, ilat, ilon2)
        @test cell_volume(g, id1) ≈ cell_volume(g, id2) atol=1e-14
    end
end

@testset "geographic comparison near equator" begin
    R = 6371.0
    g = LatLonGrid(lat_edges=[-10.0, 10.0], lon_edges=[0.0, 10.0], R=R)

    # Geographic formula: R² * (sin(lat2) - sin(lat1)) * Δlon_rad
    lat1, lat2 = deg2rad(-10.0), deg2rad(10.0)
    Δlon = deg2rad(10.0)
    geographic_area = R^2 * (sin(lat2) - sin(lat1)) * Δlon

    geodesic_area = cell_volume(g, 1)

    # Near equator, geodesic (great-circle) and geographic should be very close
    @test geodesic_area ≈ geographic_area rtol=0.01
end

@testset "polar cell volume correctness" begin
    R = 1.0
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 120.0, 240.0, 360.0], R=R)

    # South polar cell (ilat=1): should be a spherical triangle with area π/3 steradians
    # (3 cells covering the southern hemisphere, each is an octant-like triangle)
    south_vol = cell_volume(g, _cell_linear_index(g, 1, 1))
    @test south_vol > 0.0
    @test !isnan(south_vol)

    # North polar cell (ilat=2): same by symmetry
    north_vol = cell_volume(g, _cell_linear_index(g, 2, 1))
    @test north_vol ≈ south_vol atol=1e-14

    # Equatorial band total: 4π - 2 * polar area
    equatorial_total = sum(ilon -> cell_volume(g, _cell_linear_index(g, 2, ilon)), 1:g.nlon)
    @test equatorial_total ≈ 4 * π * R^2 - 2 * south_vol atol=1e-12
end
```

Note: Tests reference `ManifoldMeshes._spherical_triangle_area` and `_cell_linear_index`. These are internal functions accessed via module qualification. If this is a problem, expose them as `ManifoldMeshes._spherical_triangle_area` by ensuring the function is defined at module scope (not nested). The current code in Task 4 defines them at module scope, so this works.

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `cell_volume` undefined.

- [ ] **Step 3: Implement cell_volume**

Add to `src/sphere/latlon.jl`:

```julia
# ── Geometry: cell_volume ──────────────────────────────────────────

function cell_volume(g::LatLonGrid, cell_id::Int)
    ilat, ilon = _cell_indices(g, cell_id)
    return g.cell_volumes[ilat, ilon]
end
```

Also add `cell_centroid` (simple cache read, needed for later tasks):

```julia
# ── Geometry: cell_centroid ────────────────────────────────────────

function cell_centroid(g::LatLonGrid, cell_id::Int)
    ilat, ilon = _cell_indices(g, cell_id)
    return g.cell_centroids[ilat, ilon]
end
```

- [ ] **Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add src/sphere/latlon.jl test/test_latlon_geometry.jl
git commit -m "feat: add cell_volume and cell_centroid with lobe test"
```

---

## Task 6: Cell Centroid Tests

**Files:**
- Modify: `test/test_latlon_geometry.jl`

- [ ] **Step 1: Write additional centroid tests**

Append to `test/test_latlon_geometry.jl`:

```julia
@testset "cell_centroid basics" begin
    R = 1.0
    g = LatLonGrid(lat_edges=[-90.0, 90.0], lon_edges=[0.0, 120.0, 240.0, 360.0], R=R)

    # Single latitude band (full sphere split into 3 longitudinal sectors)
    # Each cell is a lune covering 120° longitude, centered at lat=0
    for ilon in 1:g.nlon
        c = cell_centroid(g, _cell_linear_index(g, 1, ilon))
        # Should be on the sphere surface
        @test abs(norm(c) - R) < 1e-10
        # Should be at equator (z ≈ 0)
        @test abs(c[3]) < 0.1  # loose tolerance — Frechet mean of 4 vertices
    end
end

@testset "cell_centroid on sphere surface" begin
    g = LatLonGrid(lat_edges=collect(-90.0:10.0:90.0),
                   lon_edges=collect(0.0:15.0:360.0), R=1.0)
    for i in 1:num_cells(g)
        c = cell_centroid(g, i)
        @test abs(norm(c) - 1.0) < 1e-10
    end
end
```

- [ ] **Step 2: Run tests**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: ALL PASS

- [ ] **Step 3: Commit**

```bash
git add test/test_latlon_geometry.jl
git commit -m "test: add cell_centroid correctness tests"
```

---

## Task 7: Connectivity (cell_nodes, cell_cells, node_cells, cell_edges)

**Files:**
- Create: `test/test_latlon_connectivity.jl`
- Modify: `src/sphere/latlon.jl`

- [ ] **Step 1: Write the failing test**

```julia
# test/test_latlon_connectivity.jl
@testset "cell_nodes consistency" begin
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 90.0, 180.0, 270.0, 360.0])

    # Cell 1 (ilat=1, ilon=1): SW corner cell
    sn = cell_nodes(g, 1)
    @test length(sn) == 4
    # SW node should be at south pole
    @test node_coordinates(g, sn[1]) ≈ SVector(0.0, 0.0, -1.0) atol=1e-12

    # All cells should return NTuple{4,Int}
    for i in 1:num_cells(g)
        @test cell_nodes(g, i) isa NTuple{4,Int}
    end
end

@testset "cell_cells periodicity" begin
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 90.0, 180.0, 270.0, 360.0])

    # Cell at (ilat=1, ilon=1): west neighbor is (1, nlon=4)
    _, _, west, _ = cell_cells(g, 1)
    @test west == _cell_linear_index(g, 1, 4)

    # Cell at (ilat=1, ilon=4): east neighbor is (1, 1)
    _, _, _, east = cell_cells(g, _cell_linear_index(g, 1, 4))
    @test east == _cell_linear_index(g, 1, 1)
end

@testset "cell_cells sentinel values" begin
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 360.0])

    # Bottom row (ilat=1): south neighbor = 0
    south, north, west, east = cell_cells(g, 1)
    @test south == 0
    @test north == 2  # only cell in row 2
    @test west == 1   # periodic (only 1 column)
    @test east == 1   # periodic

    # Top row (ilat=2): north neighbor = 0
    south, north, west, east = cell_cells(g, 2)
    @test south == 1
    @test north == 0

    # Always NTuple{4,Int}
    for i in 1:num_cells(g)
        @test cell_cells(g, i) isa NTuple{4,Int}
    end
end

@testset "node_cells interior vs pole" begin
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 90.0, 180.0, 270.0, 360.0])

    # Interior node (ilat=2, ilon=1): equator, lon=0 — 4 adjacent cells
    nc = node_cells(g, _node_linear_index(g, 2, 1))
    @test length(nc) == 4

    # South pole node (ilat=1, ilon=1): 2 adjacent cells
    nc = node_cells(g, _node_linear_index(g, 1, 1))
    @test length(nc) == 2

    # North pole node (ilat=3, ilon=1): 2 adjacent cells
    nc = node_cells(g, _node_linear_index(g, 3, 1))
    @test length(nc) == 2
end

@testset "cell_edges consistency with cell_nodes" begin
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 90.0, 180.0, 270.0, 360.0])

    for cell_id in 1:num_cells(g)
        edges = cell_edges(g, cell_id)
        @test edges isa NTuple{4,Int}

        # South edge endpoints should match SW and SE nodes
        sn = cell_nodes(g, cell_id)
        n1, n2 = ManifoldMeshes._edge_endpoints(g, edges[1])  # south edge
        @test n1 == node_coordinates(g, sn[1])  # SW
        @test n2 == node_coordinates(g, sn[2])  # SE
    end
end

@testset "cell_edges periodicity" begin
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 90.0, 180.0, 270.0, 360.0])

    # Cell (1,1) and cell (1,4) share the western edge of cell (1,1) / eastern edge of cell (1,4)
    _, _, w1, _ = cell_edges(g, _cell_linear_index(g, 1, 1))
    _, _, _, e4 = cell_edges(g, _cell_linear_index(g, 1, 4))
    @test w1 == e4  # same edge, periodic wrap
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `cell_nodes`, `cell_cells`, `node_cells`, `cell_edges`, `_edge_endpoints` undefined.

- [ ] **Step 3: Implement connectivity functions**

Add to `src/sphere/latlon.jl`:

```julia
# ── Connectivity ───────────────────────────────────────────────────

function cell_nodes(g::LatLonGrid, cell_id::Int)
    ilat, ilon = _cell_indices(g, cell_id)
    ilon_next = ilon == g.nlon ? 1 : ilon + 1
    return (
        _node_linear_index(g, ilat,   ilon),       # SW
        _node_linear_index(g, ilat,   ilon_next),   # SE
        _node_linear_index(g, ilat+1, ilon_next),   # NE
        _node_linear_index(g, ilat+1, ilon),         # NW
    )
end

function cell_cells(g::LatLonGrid, cell_id::Int)
    ilat, ilon = _cell_indices(g, cell_id)
    south = ilat > 1     ? _cell_linear_index(g, ilat - 1, ilon) : 0
    north = ilat < g.nlat ? _cell_linear_index(g, ilat + 1, ilon) : 0
    west  = _cell_linear_index(g, ilat, ilon == 1    ? g.nlon : ilon - 1)
    east  = _cell_linear_index(g, ilat, ilon == g.nlon ? 1      : ilon + 1)
    return (south, north, west, east)
end

function node_cells(g::LatLonGrid, node_id::Int)
    ilat, ilon = _node_indices(g, node_id)
    cells = Int[]
    for i in (ilat - 1, ilat)
        i < 1 || i > g.nlat && continue
        for j in (ilon - 1, ilon)
            jj = j < 1 ? g.nlon : (j > g.nlon ? 1 : j)
            if 1 ≤ i ≤ g.nlat && 1 ≤ jj ≤ g.nlon
                push!(cells, _cell_linear_index(g, i, jj))
            end
        end
    end
    return cells
end

function cell_edges(g::LatLonGrid, cell_id::Int)
    ilat, ilon = _cell_indices(g, cell_id)
    n_h = (g.nlat + 1) * g.nlon  # number of horizontal edges
    south = (ilat - 1) * g.nlon + ilon
    north = ilat * g.nlon + ilon
    west  = n_h + (ilat - 1) * g.nlon + ilon
    east_ilon = ilon == g.nlon ? 1 : ilon + 1
    east  = n_h + (ilat - 1) * g.nlon + east_ilon
    return (south, north, west, east)
end
```

Also add `_edge_endpoints` (needed by tests and later tasks):

```julia
# ── Internal: Edge Endpoints ───────────────────────────────────────

function _edge_endpoints(g::LatLonGrid, edge_id::Int)
    n_h = (g.nlat + 1) * g.nlon
    if edge_id <= n_h
        # Horizontal edge: along a latitude circle
        idx = edge_id - 1
        ilat = div(idx, g.nlon) + 1
        ilon = rem(idx, g.nlon) + 1
        ilon_next = ilon == g.nlon ? 1 : ilon + 1
        return (g.nodes[ilat, ilon], g.nodes[ilat, ilon_next])
    else
        # Vertical edge: along a longitude line
        idx = edge_id - n_h - 1
        ilat = div(idx, g.nlon) + 1
        ilon = rem(idx, g.nlon) + 1
        return (g.nodes[ilat, ilon], g.nodes[ilat + 1, ilon])
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add src/sphere/latlon.jl test/test_latlon_connectivity.jl
git commit -m "feat: add connectivity (cell_nodes, cell_cells, node_cells, cell_edges)"
```

---

## Task 8: Edge Operations + Boundary Markers

**Files:**
- Create: `test/test_latlon_normals.jl`
- Modify: `src/sphere/latlon.jl`

- [ ] **Step 1: Write the failing test**

```julia
# test/test_latlon_normals.jl
@testset "edge_length equatorial" begin
    R = 1.0
    g = LatLonGrid(lat_edges=[-90.0, 90.0], lon_edges=[0.0, 90.0, 180.0, 270.0, 360.0], R=R)

    # Vertical edge at equator, spanning 180° latitude: length should be πR
    # Edge connecting equator to north pole = half great circle = πR
    # Cell 1 (ilat=1, ilon=1): west edge (vertical, ilon=1)
    # v(1,1) connects node(1,1) to node(2,1) = south pole to equator
    n_h = (g.nlat + 1) * g.nlon
    west_edge = n_h + (1-1)*g.nlon + 1  # = 4 + 1 = 5
    @test edge_length(g, west_edge) ≈ π/2 * R atol=1e-12  # south pole to equator = π/2

    # Horizontal edge at equator, spanning 90° longitude
    # h(2,1) connects node(2,1) to node(2,2) = equator at lon=0 to equator at lon=90
    eq_edge = (2-1)*g.nlon + 1  # = 4 + 1 = 5... wait
    # Actually h(ilat=2, ilon=1) = (2-1)*4 + 1 = 5
    eq_edge = (2-1)*g.nlon + 1
    @test edge_length(g, eq_edge) ≈ π/2 * R atol=1e-12  # 90° on equator
end

@testset "edge_midpoint equidistant" begin
    g = LatLonGrid(lat_edges=[-90.0, 90.0], lon_edges=[0.0, 180.0, 360.0], R=1.0)

    for edge_id in 1:num_edges(g)
        mp = edge_midpoint(g, edge_id)
        n1, n2 = ManifoldMeshes._edge_endpoints(g, edge_id)
        # Midpoint should be equidistant from both endpoints
        d1 = Manifolds.distance(g.manifold, mp, n1)
        d2 = Manifolds.distance(g.manifold, mp, n2)
        @test d1 ≈ d2 atol=1e-12
        # Midpoint should be on the sphere
        @test abs(norm(mp) - 1.0) < 1e-10
    end
end

@testset "edge_outward_normal orthogonality" begin
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 90.0, 180.0, 270.0, 360.0])

    for cell_id in 1:num_cells(g)
        edges = cell_edges(g, cell_id)
        cc = cell_centroid(g, cell_id)
        for edge_id in edges
            result = edge_outward_normal(g, edge_id, cell_id)
            bp, n = result.base_point, result.normal

            # Normal should be in tangent space: dot(bp, n) ≈ 0
            @test abs(dot(bp, n)) < 1e-10

            # Normal should be unit length
            @test abs(norm(n) - 1.0) < 1e-10

            # Normal should point outward: dot(normal, centroid - base_point) > 0
            # (centroid is inside the cell, outward normal points away from cell)
            diff = cc - bp
            proj_diff = diff - dot(diff, bp) * bp  # project to tangent space
            @test dot(n, proj_diff) < 0  # normal points away from cell interior
        end
    end
end

@testset "polar collapse: degenerate edge normal" begin
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 90.0, 180.0, 270.0, 360.0])

    # South pole horizontal edge (ilat=1): all nodes at south pole, zero-length edge
    pole_edge = 1  # h(1,1) = (1-1)*4 + 1 = 1
    result = edge_outward_normal(g, pole_edge, 1)
    @test result.normal == zero(SVector{3,Float64})
    @test !any(isnan, result.normal)

    # North pole horizontal edge (ilat=nlat+1=3)
    n_h = (g.nlat + 1) * g.nlon
    north_pole_edge = n_h  # last horizontal edge
    # This edge belongs to cells in row nlat=2
    result = edge_outward_normal(g, north_pole_edge, _cell_linear_index(g, 2, 1))
    @test result.normal == zero(SVector{3,Float64})
    @test !any(isnan, result.normal)
end

@testset "boundary markers (full sphere = no boundary)" begin
    g = LatLonGrid(lat_edges=[-90.0, 90.0], lon_edges=[0.0, 360.0])
    @test boundary_nodes(g, :any) == Int[]
    @test boundary_edges(g, :any) == Int[]
    @test boundary_nodes(g, 0) == Int[]
    @test boundary_edges(g, 0) == Int[]
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `edge_length`, `edge_midpoint`, `edge_outward_normal`, `boundary_nodes`, `boundary_edges` undefined.

- [ ] **Step 3: Implement edge operations and boundary markers**

Add to `src/sphere/latlon.jl`:

```julia
# ── Geometry: edge_length ──────────────────────────────────────────

function edge_length(g::LatLonGrid, edge_id::Int)
    n1, n2 = _edge_endpoints(g, edge_id)
    return Manifolds.distance(g.manifold, n1, n2)
end

# ── Geometry: edge_midpoint ───────────────────────────────────────

function edge_midpoint(g::LatLonGrid, edge_id::Int)
    n1, n2 = _edge_endpoints(g, edge_id)
    if Manifolds.distance(g.manifold, n1, n2) < 1e-14
        return SVector{3,Float64}(n1)
    end
    return SVector{3,Float64}(Manifolds.mid_point(g.manifold, n1, n2))
end

# ── Geometry: edge_outward_normal ──────────────────────────────────

function edge_outward_normal(g::LatLonGrid, edge_id::Int, cell_id::Int)
    n1, n2 = _edge_endpoints(g, edge_id)

    # Guard: polar collapse — degenerate edge with coincident endpoints
    if Manifolds.distance(g.manifold, n1, n2) < 1e-14
        midpoint = n1
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

# ── Boundary Markers ──────────────────────────────────────────────

boundary_nodes(g::LatLonGrid, marker) = Int[]
boundary_edges(g::LatLonGrid, marker) = Int[]
```

- [ ] **Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add src/sphere/latlon.jl test/test_latlon_normals.jl
git commit -m "feat: add edge operations (length, midpoint, outward normal) and boundary markers"
```

---

## Self-Review Checklist

### Spec Coverage

| Spec Section | Task |
|---|---|
| 1. Design Philosophy | Architectural principle followed throughout |
| 2. Design Decisions | All choices reflected in code (SVector, geographic lat, pre-computation) |
| 3. Dependencies | Project.toml + `using` in entry point |
| 4.1 TopologyStyle | Task 2 |
| 4.2 AbstractLocation | Task 2 |
| 4.3 AbstractManifoldMesh | Task 3 |
| 5.1 Global Info (manifold, num_cells/nodes/edges) | Task 4 |
| 5.2 Connectivity (cell_nodes, cell_cells, cell_edges, node_cells) | Task 7 |
| 5.3 Geometry (node_coordinates, cell_centroid, cell_volume, edge_length, edge_midpoint, edge_outward_normal) | Tasks 4, 5, 6, 8 |
| 5.4 Boundary Markers | Task 8 |
| 6.1 Struct Definition | Task 4 |
| 6.2 Constructor | Task 4 |
| 6.3 Cell Indexing | Task 4 (index helpers) |
| 6.4 Coordinate Conversion | Task 4 (node_coordinates) |
| 6.5 Cell Volume | Task 5 (_spherical_triangle_area) |
| 6.6 Edge Length | Task 8 |
| 6.7 Edge Outward Normal | Task 8 (with polar collapse guard) |
| 6.8 Connectivity | Task 7 (with sentinel values, periodicity) |
| 6.9 Boundary Markers | Task 8 |
| 6.10 TopologyStyle | Task 4 |
| 8.1 Lobe Test | Task 5 |
| 8.2 Symmetry Test | Task 5 |
| 8.3 Connectivity Periodicity | Task 7 |
| 8.4 Normal Orthogonality | Task 8 |
| 8.5 Volume vs Geographic | Task 5 |
| 8.6 Polar Collapse | Task 8 |

**Coverage: complete.** All spec sections map to tasks.

### Placeholder Scan

No TBD, TODO, "similar to", or codeless steps found. Every step has complete code, exact file paths, and exact commands.

### Type Consistency

- `_cell_indices`, `_node_indices`, `_cell_linear_index`, `_node_linear_index` — consistent naming throughout
- `_spherical_triangle_area` — defined once in Task 4, referenced in Task 5 tests
- `_edge_endpoints` — defined in Task 7, referenced in Task 8 tests
- `cell_cells` returns `NTuple{4,Int}` with `0` sentinel — consistent with spec
- `edge_outward_normal` returns `NamedTuple{(:base_point, :normal), ...}` — consistent with spec
- All coordinate returns are `SVector{3,Float64}` — consistent with spec
