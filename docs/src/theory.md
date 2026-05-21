# Theory Background

This page explains the mathematics behind ManifoldMeshes.jl, progressing from intuition to rigorous detail.

## Part 1: Intuitive Explanation

### What Is a Manifold Mesh?

Imagine covering the surface of a ball with tiles. Each tile is a small region bounded by curves on the surface, and the tiles fit together edge-to-edge with no gaps or overlaps. A **manifold mesh** is exactly this: a partition of a smooth curved surface (the manifold) into discrete cells connected by shared edges and nodes.

In ManifoldMeshes.jl, the manifold is the unit sphere $S^2$, and the tiles are spherical quadrilaterals whose boundaries are **geodesic arcs** -- segments of great circles, the straightest possible curves on a sphere.

### Why Geodesic Cells?

A natural first idea is to carve the sphere along lines of constant latitude and longitude. But coordinate lines of latitude converge at the poles: cells near the equator are wide, while cells near the poles degenerate into thin wedges. This **polar singularity** distorts cell shapes and causes numerical problems in finite-volume methods.

Geodesic arcs avoid this. A great-circle arc between two points on the sphere is the unique shortest path connecting them -- the generalization of a straight line to curved geometry. When cell boundaries are geodesic arcs, every cell has well-defined, well-behaved edges regardless of its position on the sphere. The cells are geometrically meaningful objects on the manifold itself, not artifacts of a particular coordinate chart.

### Why Separate Geometry from Topology?

The **topology** of a mesh -- which cells are neighbors, which nodes bound which cells, how cells connect -- is purely combinatorial. The **geometry** -- where nodes sit on the manifold, how long edges are, how large cells are -- depends on the manifold's metric.

By keeping these concerns separate, the same topological structure (e.g., a regular grid) can be realized on different manifolds without changing the connectivity code. All continuous mathematics (distances, geodesics, projections) is delegated to Manifolds.jl. This library never re-implements a geodesic.

## Part 2: Grid Construction Principles

ManifoldMeshes.jl provides four grid types on $S^2$, each constructed by a different strategy. All four share the same interface (see the [API Reference](api.md)) and produce cells whose boundaries are geodesic arcs.

### LatLonGrid

The most direct construction: subdivide the colatitude $\theta \in [0, \pi]$ into `nlat` bands and the longitude $\phi \in [0, 2\pi)$ into `nlon` sectors. The intersection of a latitude band and a longitude sector defines a cell.

Each cell is a spherical quadrilateral. Its four vertices are the $(\theta, \phi)$ grid points, embedded as Cartesian points on $S^2$:

$$x = (\sin\theta\cos\phi,\; \sin\theta\sin\phi,\; \cos\theta).$$

**Key property:** Structured Cartesian indexing. Cell `(i, j)` has neighbors `(i-1, j)`, `(i+1, j)`, `(i, j-1)`, `(i, j+1)` with periodic wrapping in longitude.

**Limitation:** Cells near the poles become extremely narrow because all longitude lines converge at $\theta = 0$ and $\theta = \pi$. The polar "cells" degenerate into triangles or even lunes. This grid is best suited for applications where the polar singularity can be tolerated or where compatibility with legacy latitude-longitude data formats is required.

### CubedSphereGrid

To eliminate polar singularities, project the six faces of an inscribed cube onto the sphere. Each cube face is a square domain $[-1, 1]^2$ with local coordinates $(s, t)$. The **gnomonic projection** maps a point on the cube face to the sphere by drawing a ray from the sphere's center through the cube point:

$$\mathbf{x}(s, t) = \frac{1}{\sqrt{1 + s^2 + t^2}}\,(v_1 s + v_2 t + v_3),$$

where $(v_1, v_2, v_3)$ are orthonormal vectors defining the face orientation. Each face is subdivided into $n \times n$ cells, giving $6n^2$ total cells.

**Key property:** No coordinate singularity. Every face has uniform quasi-regular cells. The trade-off is mild area variation across faces and the need to handle cross-face topology explicitly. The `MultiPatch` trait and [`cell_face`](@ref) / [`cell_local_2d`](@ref) functions manage this multi-patch structure.

### ReducedGaussianGrid

Used in global spectral weather models (e.g., IFS, GFS). Latitude bands are placed at the **roots of Legendre polynomials** -- the Gaussian quadrature points -- which are optimal for spectral transforms (Legendre expansions in latitude, Fourier expansions in longitude). The number of longitude points per band decreases toward the poles following an octahedral reduction pattern, producing cells of approximately equal area.

**Key property:** Optimality for spectral methods. The Gaussian latitudes are the nodes of Gaussian quadrature, making trapezoidal integration in latitude exact for polynomials up to a known degree. The variable longitude count avoids the over-sampling at high latitudes that plagues regular LatLonGrids.

### HEALPixGrid

The **H**ierarchical **E**qual **A**rea iso-**L**atitude **Pix**elization divides the sphere into 12 base pixels (4 equatorial diamonds and 8 polar triangles), each recursively subdivided into $n_{\mathrm{side}} \times n_{\mathrm{side}}$ sub-pixels, giving $12\,n_{\mathrm{side}}^2$ total cells.

**Key property:** Exact equal-area property. Every cell has area $4\pi R^2 / (12\,n_{\mathrm{side}}^2)$ regardless of position on the sphere. This makes HEALPix ideal for all-sky statistical analyses, pixel-based CMB mapping, and any application where uniform solid-angle sampling is required. Cells are also arranged on iso-latitude rings, enabling fast spherical harmonic transforms via the `:ring` ordering scheme.

## Part 3: Mathematical Foundations

### The Sphere as a Riemannian Manifold

The unit sphere $S^2 = \{\mathbf{x} \in \mathbb{R}^3 : \|\mathbf{x}\| = 1\}$ is a two-dimensional Riemannian manifold with the round metric inherited from the Euclidean embedding. In spherical coordinates $(\theta, \phi)$ where $\theta$ is colatitude and $\phi$ is longitude, the metric tensor is

$$g = \mathrm{d}\theta \otimes \mathrm{d}\theta + \sin^2\!\theta\;\mathrm{d}\phi \otimes \mathrm{d}\phi,$$

and the Riemannian volume form is

$$\mathrm{d}A = \sin\theta\;\mathrm{d}\theta \wedge \mathrm{d}\phi.$$

The total surface area is $\int_{S^2} \mathrm{d}A = 4\pi$.

### Geodesics and Cell Boundaries

A **geodesic** on $S^2$ is a great circle: the intersection of $S^2$ with a plane through the origin. The unique shortest geodesic from point $\mathbf{a}$ to point $\mathbf{b}$ traces the arc of a great circle in the plane spanned by $\mathbf{a}$ and $\mathbf{b}$.

In ManifoldMeshes.jl, every cell boundary edge is a geodesic arc. The geodesic arc length between two nodes is computed via Manifolds.jl:

$$d(\mathbf{a}, \mathbf{b}) = R \cdot \arccos\!\left(\frac{\mathbf{a} \cdot \mathbf{b}}{R^2}\right),$$

which [`edge_length`](@ref) returns directly.

### Cell Volume: Spherical Excess via l'Huilier's Formula

Each quadrilateral cell is split into two spherical triangles along a diagonal, and the areas of the two triangles are summed. The function `_spherical_triangle_area` uses **l'Huilier's formula**, which expresses the area of a spherical triangle in terms of its side lengths.

Given a spherical triangle with side lengths $a$, $b$, $c$ (in radians on a unit sphere) and semi-perimeter $s = (a + b + c) / 2$, the **spherical excess** $E$ is

$$\tan\frac{E}{4} = \sqrt{\tan\frac{s}{2}\,\tan\frac{s-a}{2}\,\tan\frac{s-b}{2}\,\tan\frac{s-c}{2}}.$$

The triangle area on a sphere of radius $R$ is then $A = R^2 E$. The side lengths are computed from the dot products of the vertex positions:

$$a = \arccos\!\left(\frac{\mathbf{B} \cdot \mathbf{C}}{R^2}\right), \quad b = \arccos\!\left(\frac{\mathbf{A} \cdot \mathbf{C}}{R^2}\right), \quad c = \arccos\!\left(\frac{\mathbf{A} \cdot \mathbf{B}}{R^2}\right).$$

For polar cells where one diagonal degenerates (e.g., the pole), the code falls back to the other diagonal, or to a lune formula:

$$A_{\mathrm{lune}} = R^2\,(\sin\theta_2 - \sin\theta_1)\,|\Delta\phi|.$$

A fundamental sanity check: $\sum_{\text{all cells}} \mathrm{cell\_volume} = 4\pi R^2$.

### Edge Outward Normal: Tangent-Space Semantics

The function [`edge_outward_normal`](@ref) returns a `NamedTuple{(:base_point, :normal)}`. The `base_point` is the geodesic midpoint of the edge (computed via `Manifolds.mid_point`), and `normal` is a unit vector in the **tangent space** at that point, perpendicular to the edge and pointing away from the reference cell.

The tangent space $T_{\mathbf{p}} S^2$ at a point $\mathbf{p}$ is the plane perpendicular to $\mathbf{p}$. A vector $\mathbf{v} \in \mathbb{R}^3$ is in the tangent space if and only if $\mathbf{p} \cdot \mathbf{v} = 0$. Manifolds.jl provides projection operators that ensure the computed normal lies correctly in this tangent plane.

This tangent-space representation is essential for finite-volume flux computations: the flux across an edge is $\mathbf{F} \cdot \hat{\mathbf{n}}$ evaluated at the edge midpoint, where $\hat{\mathbf{n}}$ is the outward normal in the tangent space.

### Relationship to Manifolds.jl

ManifoldMeshes.jl delegates all continuous differential geometry to [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/). Specifically:

| ManifoldMeshes function | Manifolds.jl operation |
|---|---|
| [`edge_length`](@ref) | `Manifolds.distance(M, p, q)` |
| [`edge_midpoint`](@ref) | `Manifolds.mid_point(M, p, q)` |
| [`cell_centroid`](@ref) | `Manifolds.mean(M, vertices)` (Frechet mean) |
| [`edge_outward_normal`](@ref) | `Manifolds.project(M, base_point, vector)` |

This library never re-implements a geodesic, never recomputes a distance formula, and never derives a tangent-space projection from scratch.

## Part 4: Design Principles

### Geometry Belongs to Manifolds, Topology Belongs to Meshes

This is the central architectural tenet. The mesh knows *which* cells are neighbors, *which* nodes bound a cell, and *how* cells are indexed. It does not know *how far apart* two points are or *what curve* connects them -- that is the manifold's job. By delegating geometry to Manifolds.jl, the mesh code stays focused on combinatorial structure and can in principle be ported to other manifolds (e.g., the torus, hyperbolic surfaces) by changing only the manifold parameter.

### No Physics

A mesh struct in ManifoldMeshes.jl never holds field data -- no temperatures, no wind vectors, no chemical concentrations. Physical quantities belong to a separate library (ManifoldFields.jl). This separation keeps the mesh layer small, testable, and reusable across applications.

### Performance Through Caching

All geometric quantities -- cell volumes, centroids, node positions -- are pre-computed at grid construction time and stored in dense arrays. Queries like [`cell_volume`](@ref) and [`cell_centroid`](@ref) are simple array lookups with no runtime computation. This trades memory for speed, which is the right trade-off for simulation codes that query these values millions of times per time step.

### Pure Manifold School

Cell boundaries are geodesic arcs (great circles on $S^2$), not coordinate lines. This means a cell's shape is determined by the intrinsic geometry of the manifold, not by the choice of coordinate chart. The advantage is geometric consistency: cell volumes, edge lengths, and normals are all computed within the same Riemannian framework.

### Trait System for Compile-Time Dispatch

The type traits -- `TopologyStyle` (`IsGrid`, `IsSemiGrid`, `IsMesh`), `CellTypeStyle` (`IsUniform{K}`, `IsMixed{MAX_K}`), `PatchStyle` (`NoPatch`, `MultiPatch{N}`), and `AbstractLocation` (`NodeLoc`, `CellLoc`, `EdgeLoc`) -- enable the Julia compiler to specialize code paths at compile time. A loop over cells of an `IsGrid` mesh can use fast Cartesian indexing, while an `IsSemiGrid` mesh uses a connectivity table. The caller writes generic code against the interface; the compiler picks the fastest path.
