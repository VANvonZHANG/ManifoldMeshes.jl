using Manifolds: Sphere
using StaticArrays: SMatrix

# -- Struct --

struct HEALPixGrid{M <: AbstractManifold} <: AbstractManifoldMesh{M}
    manifold::M
    nside::Int
    R::Float64
    ordering::Symbol
    nodes::Vector{SVector{3, Float64}}
    cell_volumes::Vector{Float64}
    cell_centroids::Vector{SVector{3, Float64}}
    _cell_nodes::Vector{NTuple{4, Int}}
    _cell_edges::Vector{NTuple{4, Int}}
    _cell_cells::Vector{NTuple{4, Int}}
    _edge_nodes::Vector{NTuple{2, Int}}
    _dual::Base.RefValue{Union{Nothing, AbstractManifoldMesh{M}}}
    rotation::SMatrix{3, 3, Float64, 9}
end

# -- Internal: Ring-to-Nested Permutation --

# Standard HEALPix base pixel layout (0-indexed):
#      0  1
#   2  3  4  5
#   6  7  8  9
#     10 11

# Ring ranges for each base pixel (inclusive, 1-indexed ring numbers)
# For nside=1, each base pixel is exactly one ring cell.
# For nside>1, each base pixel spans multiple rings.
# North polar cap base pixels (0,1,2,3) span rings 1:nside
# Equatorial base pixels (4,5,6,7,8,9,10,11) span rings (nside+1):3nside
# South polar cap base pixels (8,9,10,11) span rings (3nside+1):(4nside-1)
# Wait, that's not right either. Let me think more carefully.
#
# Actually, the base pixels are arranged in a specific pattern on the sphere.
# For ring ordering, the cells are ordered by ring, then by longitude within each ring.
# For nested ordering, the cells are ordered by base pixel, then by Morton index within each base pixel.
#
# The key insight is: for a given nside, we can compute the (theta, phi) of each nested cell
# using the standard HEALPix formulas, then find the closest ring cell.
#
# Standard HEALPix nested cell -> (theta, phi) formulas:
# - Base pixel determines the coarse region
# - Within a base pixel, (ix, iy) in [0, nside-1] determines the fine position
# - The nested index is: base * nside^2 + morton_index(ix, iy)
#   where base is 0-indexed and morton_index is also 0-indexed.

"""
    _morton_decode(m::Int) -> (ix::Int, iy::Int)

Decode a Morton (Z-order) index into 2D grid coordinates.
Interleaves bits: m = ...b2a2b1a1b0a0 -> ix = ...a2a1a0, iy = ...b2b1b0.
"""
function _morton_decode(m::Int)
    ix = 0
    iy = 0
    bit = 0
    while m > 0
        ix |= (m & 1) << bit
        m >>= 1
        iy |= (m & 1) << bit
        m >>= 1
        bit += 1
    end
    return (ix, iy)
end

"""
    _morton_encode(ix::Int, iy::Int) -> m::Int

Encode 2D grid coordinates into a Morton (Z-order) index.
"""
function _morton_encode(ix::Int, iy::Int)
    m = 0
    bit = 0
    while ix > 0 || iy > 0
        m |= (ix & 1) << (2 * bit)
        m |= (iy & 1) << (2 * bit + 1)
        ix >>= 1
        iy >>= 1
        bit += 1
    end
    return m
end

"""
    _nested_cell_center(nside, base, ix, iy) -> SVector{3,Float64}

Compute the 3D Cartesian coordinates of a HEALPix cell center in nested ordering.

Parameters:
- nside: HEALPix resolution parameter
- base: base pixel index (0-11)
- ix, iy: integer coordinates within the base pixel, in [0, nside-1]

Returns the cell center as a unit vector (R=1).

This implements the standard HEALPix nested ordering geometry.
"""
function _nested_cell_center(nside::Int, base::Int, ix::Int, iy::Int)
    # Normalize coordinates to [0, 1] within the base pixel
    # The HEALPix cell centers are at (ix + 0.5) / nside
    u = (ix + 0.5) / nside
    v = (iy + 0.5) / nside

    # Base pixel layout on the sphere (standard HEALPix):
    # The 12 base pixels correspond to the faces of a rhombic dodecahedron.
    # We use a simplified approach: for each base pixel, define a local coordinate
    # system and map (u, v) to spherical coordinates.
    #
    # The standard approach is to use the HEALPix projection (H=1.5, X=1.0):
    # For a given base pixel, compute (x, y) in the HEALPix projection plane,
    # then project back to the sphere.
    #
    # Simplified approach for our purposes:
    # We'll compute the nested cell center directly from the standard HEALPix formulas
    # by determining which ring the cell belongs to and its position within that ring.

    # For the nested ordering, cells within a base pixel are arranged in a grid.
    # The base pixels have different orientations on the sphere.
    #
    # Let's use the direct approach: compute the nested index, then use the
    # standard HEALPix formulas to get (theta, phi).

    # Convert (base, ix, iy) to nested index
    morton = _morton_encode(ix, iy)
    nested_idx = base * nside * nside + morton  # 0-indexed

    # Now convert nested index to (theta, phi) using standard HEALPix formulas
    return _nested_to_ang(nside, nested_idx)
end

"""
    _nested_to_ang(nside, nested_idx) -> SVector{3,Float64}

Convert a 0-indexed nested cell index to 3D Cartesian coordinates on the unit sphere.
Implements the standard HEALPix nested ordering geometry.
"""
function _nested_to_ang(nside::Int, nested_idx::Int)
    n_cells = 12 * nside * nside
    @assert 0 <= nested_idx < n_cells

    npface = nside * nside
    face = div(nested_idx, npface)  # base pixel (0-11)
    ipf = mod(nested_idx, npface)   # index within face

    # Decode Morton index to (ix, iy) within the face
    ix, iy = _morton_decode(ipf)

    # Standard HEALPix formulas for nested ordering
    # From the HEALPix paper (Gorski et al. 2005)
    #
    # For each face, we compute (x, y) in the HEALPix projection plane:
    # x = (ix + 0.5) / nside
    # y = (iy + 0.5) / nside
    #
    # Then map to (theta, phi) based on the face.

    # JR = ix + iy  # ring offset within face (0 to 2*nside-2)
    # JR = 2*nside - 2 - JR for southward faces

    # Let's use the standard implementation approach.
    # The 12 faces are arranged as:
    #   0  1
    # 2 3 4 5
    # 6 7 8 9
    #   10 11
    #
    # Faces 0-3: north polar cap
    # Faces 4-7: equatorial region
    # Faces 8-11: south polar cap

    # Actually, the standard layout is:
    # North polar cap: faces 0, 1, 2, 3
    # Equatorial: faces 4, 5, 6, 7, 8, 9, 10, 11
    # South polar cap: faces 8, 9, 10, 11... wait, that's overlapping.
    #
    # Correct standard layout:
    #   0  1
    # 2  3  4  5
    # 6  7  8  9
    #   10 11
    #
    # North cap: 0, 1, 2, 3
    # Equator: 4, 5, 6, 7, 8, 9, 10, 11
    # South cap: 8, 9, 10, 11... no, that's wrong.
    #
    # Let me look at the standard more carefully:
    # The 12 base pixels are the faces of a rhombic dodecahedron.
    # In the standard HEALPix paper, the faces are numbered:
    # North polar: 0, 1, 2, 3
    # Equatorial: 4, 5, 6, 7, 8, 9, 10, 11
    # South polar: 8, 9, 10, 11... no, the equatorial faces wrap around.
    #
    # Actually, looking at the standard implementation:
    # face 0-3: north polar region
    # face 4-7: equatorial region, northern row
    # face 8-11: equatorial region, southern row + south polar
    #
    # Hmm, let me just use the standard formulas from the HEALPix C++ code.
    # The key is that for nested ordering, we can compute (theta, phi) as follows:

    jr = ix + iy + 1  # 1-indexed ring within face (1 to 2*nside)
    # Actually, let's use the standard formulas more carefully.

    # From the HEALPix C++ implementation (healpix_base.cc):
    # For nested ordering:
    #   int nl4 = 4 * nside;
    #   int jr = (jrll[face] * nside) - ix - iy - 1;
    #   ...
    # where jrll is an array giving the ring offset for each face.

    # Standard jrll values for the 12 faces:
    jrll = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]
    jpll = [1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7]

    # Wait, that's for nside=1. For general nside, the formulas are more complex.
    # Let me use the standard approach from the HEALPix reference implementation.

    # From healpix_base.cc in the official HEALPix C++ code:
    # void Healpix_Base::nest2ring(int nside, int pix, int &ipring)
    # {
    #   int npface = nside * nside;
    #   int face = pix / npface;
    #   int ipf = pix % npface;
    #   int ix, iy;
    #   morton_decode(ipf, ix, iy);
    #   int jr = (jrll[face] * nside) - ix - iy - 1;
    #   ...
    # }
    #
    # Where jrll = {2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4}
    # and jpll = {1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7}

    # For our purposes, we want to compute (theta, phi) from (face, ix, iy).
    # The ring number (1-indexed, from north pole) is:
    # jr = jrll[face] * nside - ix - iy
    # Wait, the C++ code has -1 because it's 0-indexed internally.

    # Let me be more careful. In the C++ code:
    # jr = (jrll[face] * nside) - ix - iy - 1;
    # where ix, iy are 0-indexed.
    # This gives a 0-indexed ring number from the north pole.
    #
    # For face in north polar cap (0-3):
    #   jrll[face] = 2
    #   jr = 2*nside - ix - iy - 1
    #   When ix=iy=0: jr = 2*nside - 1 (southmost ring in face)
    #   When ix=iy=nside-1: jr = 2*nside - 2*(nside-1) - 1 = 1 (northmost ring, near pole)
    #
    # Hmm, that seems backwards. Let me check the C++ code more carefully.
    #
    # Actually, looking at the HEALPix C++ code:
    # static const int jrll[] = {2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
    # static const int jpll[] = {1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7};
    #
    # In nest2ring:
    # int jr = (jrll[face] * nside) - ix - iy - 1;
    #
    # For north polar faces (0-3), jrll = 2:
    # jr ranges from 2*nside - 2*(nside-1) - 1 = 1 to 2*nside - 1
    # Wait, that's 1 to 2*nside-1, which covers the north polar cap.
    # The north polar cap has rings 1 to nside, so this doesn't match.
    #
    # Oh wait, I think I'm confusing the indexing. Let me re-read the C++ code.
    # In the C++ code, jr is the "ring number" but it's measured from the
    # equator or something. Let me look at the full nest2ring function.
    #
    # From healpix_base.cc:
    # void Healpix_Base::nest2ring(int nside, int pix, int &ipring)
    # {
    #   int npface = nside * nside;
    #   int face = pix / npface;
    #   int ipf = pix % npface;
    #   int ix, iy;
    #   morton_decode(ipf, ix, iy);
    #   int jr = (jrll[face] * nside) - ix - iy - 1;
    #
    #   int nr, kshift, n_before;
    #   if (jr < nside)
    #   {
    #     nr = jr;
    #     n_before = 2 * nr * (nr - 1);
    #     kshift = 0;
    #   }
    #   else if (jr > 3 * nside)
    #   {
    #     nr = 4 * nside - jr;
    #     n_before = 12 * nside * nside - 2 * (nr + 1) * nr;
    #     kshift = 0;
    #   }
    #   else
    #   {
    #     nr = nside;
    #     n_before = 2 * nside * (nside - 1) + (jr - nside) * 4 * nside;
    #     kshift = (jr - nside) & 1;
    #   }
    #
    #   int jp = (jpll[face] * nr + ix - iy + 1 + kshift) / 2;
    #   if (jp > nl4) jp -= nl4;
    #   if (jp < 1) jp += nl4;
    #
    #   ipring = n_before + jp - 1;
    # }
    #
    # OK so jr is a pseudo-ring number. For north polar cap (jr < nside),
    # the actual ring number is jr (1-indexed from north pole).
    # For equatorial (nside <= jr <= 3*nside), the ring number is jr.
    # For south polar cap (jr > 3*nside), the ring number is 4*nside - jr.
    #
    # So for our purposes:
    # - If jr < nside: ring = jr, nr = jr (cells in this ring)
    # - If nside <= jr <= 3*nside: ring = jr, nr = nside
    # - If jr > 3*nside: ring = 4*nside - jr, nr = 4*nside - jr
    #
    # Wait no, looking more carefully:
    # For jr < nside: nr = jr, and this is the number of cells in the ring.
    # But the actual ring number from the north pole is... let me think.
    #
    # The total cells in the north polar cap (rings 1 to nside-1) is:
    # sum_{r=1}^{nside-1} 4*r = 2*(nside-1)*nside
    #
    # For jr < nside: n_before = 2*nr*(nr-1) where nr = jr
    # This is the number of cells before ring jr in the north polar cap.
    # So ring jr has 4*jr cells, and the ring number is jr.
    #
    # For nside <= jr <= 3*nside:
    # n_before = 2*nside*(nside-1) + (jr - nside)*4*nside
    # The first term is the total cells in the north polar cap.
    # The second term is (jr - nside) * 4 * nside, which is the number of
    # cells in the equatorial region before ring jr.
    # So ring jr has 4*nside cells, and the ring number is jr.
    #
    # For jr > 3*nside:
    # nr = 4*nside - jr
    # n_before = 12*nside*nside - 2*(nr+1)*nr
    # This is the number of cells before the south polar cap ring.
    # The ring number from the south pole is nr, so from the north pole it's
    # 4*nside - 1 - (nr - 1) = 4*nside - nr = jr.
    # Wait, that gives ring = jr again. But jr > 3*nside and the total rings
    # is 4*nside - 1, so jr ranges from 3*nside+1 to 4*nside-1.
    #
    # Hmm, but the south polar cap rings are numbered from the south pole.
    # Ring 4*nside-1 is the ring just above the south pole (1 cell? No, 4 cells).
    # Actually, the south polar cap has rings numbered from 3*nside+1 to 4*nside-1,
    # with 4*(4*nside - jr) cells in ring jr.
    #
    # Let me just use the standard formulas directly.

    jrll_arr = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]
    jpll_arr = [1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7]

    jr = jrll_arr[face + 1] * nside - ix - iy - 1  # 0-indexed pseudo-ring
    # Note: in C++ this is 0-indexed, but the formulas use it as if it's 1-indexed
    # Let me add 1 to make it 1-indexed like the C++ code expects.
    # Actually no, the C++ code uses 0-indexed jr in the comparisons.
    # Let me re-check: in C++, jr = (jrll[face] * nside) - ix - iy - 1
    # For face=0, nside=1, ix=0, iy=0: jr = 2*1 - 0 - 0 - 1 = 1
    # Hmm, that's 1-indexed. But then the comparison is jr < nside, which for nside=1
    # would be 1 < 1, false. So it goes to the else branch.
    #
    # Wait, but for nside=1, the north polar cap has only ring 1, which has 4 cells.
    # And face 0 should be in the north polar cap. So jr < nside should be true.
    # This means jr should be 0-indexed: jr = 2*1 - 0 - 0 - 1 = 1, but we want 0.
    #
    # Hmm, maybe the C++ code uses 1-indexed jr and the condition is jr <= nside?
    # Let me check the actual C++ code again.
    #
    # Oh I see, in the C++ code, the condition is:
    # if (jr < nside)
    # For nside=1, face=0, ix=0, iy=0: jr = 2 - 0 - 0 - 1 = 1
    # 1 < 1 is false, so it goes to else if (jr > 3*nside), which is also false.
    # Then else: nr = nside, n_before = 2*nside*(nside-1) + (jr-nside)*4*nside
    # = 0 + 0 = 0. So the first cell in face 0 is ring cell 0.
    #
    # But for nside=1, ring 1 has 4 cells, and the north polar cap has 4 cells total.
    # So face 0 cell 0 should be ring cell 0, which matches.
    #
    # OK so jr is 1-indexed in the C++ code. Let me just use the same formulas.

    # jr is already 1-indexed from the formula above
    nl4 = 4 * nside

    if jr < nside
        nr = jr
        n_before = 2 * nr * (nr - 1)
        kshift = 0
    elseif jr > 3 * nside
        nr = nl4 - jr
        n_before = n_cells - 2 * (nr + 1) * nr
        kshift = 0
    else
        nr = nside
        n_before = 2 * nside * (nside - 1) + (jr - nside) * nl4
        kshift = mod(jr - nside, 2)
    end

    jp = div(jpll_arr[face + 1] * nr + ix - iy + 1 + kshift, 2)
    if jp > nl4
        jp -= nl4
    end
    if jp < 1
        jp += nl4
    end

    # ipring is 0-indexed ring cell ID
    ipring = n_before + jp - 1

    # Now convert ipring (0-indexed) to (theta, phi)
    return _ring_to_ang(nside, ipring)
end

"""
    _ring_to_ang(nside, ipring) -> SVector{3,Float64}

Convert a 0-indexed ring cell index to 3D Cartesian coordinates on the unit sphere.
"""
function _ring_to_ang(nside::Int, ipring::Int)
    n_cells = 12 * nside * nside
    @assert 0 <= ipring < n_cells

    # Find the ring
    # North polar cap: rings 1 to nside, ring r has 4*r cells
    # Equatorial: rings nside+1 to 3*nside, each has 4*nside cells
    # South polar cap: rings 3*nside+1 to 4*nside-1, ring r has 4*(4*nside-r) cells

    n_rings = 4 * nside - 1

    # Find ring by accumulating
    ring = 1
    remaining = ipring
    while ring <= n_rings
        if ring <= nside
            n_in_ring = 4 * ring
        elseif ring <= 3 * nside
            n_in_ring = 4 * nside
        else
            n_in_ring = 4 * (4 * nside - ring)
        end

        if remaining < n_in_ring
            break
        end
        remaining -= n_in_ring
        ring += 1
    end

    # remaining is the 0-indexed cell within the ring
    i_in_ring = remaining

    # Compute theta (colatitude)
    if ring <= nside
        # North polar cap
        cos_theta = 1.0 - ring^2 / (3.0 * nside^2)
    elseif ring <= 3 * nside
        # Equatorial belt
        cos_theta = (4.0 * nside - 2.0 * ring) / (3.0 * nside)
    else
        # South polar cap
        j = n_rings - ring + 1
        cos_theta = -(1.0 - j^2 / (3.0 * nside^2))
    end
    theta = acos(cos_theta)

    # Compute phi (longitude)
    if ring <= nside
        n_in_ring = 4 * ring
    elseif ring <= 3 * nside
        n_in_ring = 4 * nside
    else
        n_in_ring = 4 * (4 * nside - ring)
    end
    phi = 2π * (i_in_ring + 0.5) / n_in_ring

    x = sin(theta) * cos(phi)
    y = sin(theta) * sin(phi)
    z = cos(theta)

    return SVector(x, y, z)
end

"""
    _ring_to_nested_permutation(nside, ring_centers_flat)

Compute the permutation mapping nested cell IDs to ring cell IDs.
Returns `perm` where `perm[nested_id] = ring_id`.

Uses the standard HEALPix nested-to-ring conversion to compute the mapping.
"""
function _ring_to_nested_permutation(nside::Int, ring_centers_flat::Vector{SVector{3, Float64}})
    n_cells = 12 * nside * nside
    perm = Vector{Int}(undef, n_cells)
    used = falses(n_cells)

    for nested_idx in 0:(n_cells - 1)
        # Compute the ring cell index for this nested cell
        p_nested = _nested_to_ang(nside, nested_idx)

        # Find the closest ring center
        best_dist = Inf
        best_ring = 0
        for (rid, rc) in enumerate(ring_centers_flat)
            d = norm(p_nested - rc)
            if d < best_dist
                best_dist = d
                best_ring = rid
            end
        end

        @assert best_ring > 0 "Failed to match nested cell $(nested_idx+1) to a ring cell"
        @assert !used[best_ring] "Ring cell $best_ring already matched by nested $(nested_idx+1)"

        perm[nested_idx + 1] = best_ring
        used[best_ring] = true
    end

    return perm
end

# -- Internal: Bounds Checking --

@inline function _check_cell_id(g::HEALPixGrid, cell_id::Int)
    @boundscheck 1 <= cell_id <= num_cells(g) ||
                 throw(BoundsError("cell_id $cell_id out of range [1, $(num_cells(g))]"))
    nothing
end

@inline function _check_node_id(g::HEALPixGrid, node_id::Int)
    @boundscheck 1 <= node_id <= num_nodes(g) ||
                 throw(BoundsError("node_id $node_id out of range [1, $(num_nodes(g))]"))
    nothing
end

# -- Constructor --

"""
    HEALPixGrid(; nside::Int, ordering::Symbol = :ring, rotation = I, R::Float64 = 1.0)

Construct a HEALPix grid on the sphere.

# Arguments
- `nside`: Resolution parameter (≥ 1). Cell count = 12 × nside²
- `ordering`: `:ring` (default) or `:nested`
- `rotation`: 3×3 rotation matrix applied to all nodes
- `R`: Sphere radius (default 1.0)

HEALPix provides an equal-area hierarchical subdivision of the sphere
into 12 base pixels, each divided into nside² quadrilateral cells.
"""
function HEALPixGrid(; nside::Int, ordering::Symbol = :ring,
        rotation = SMatrix{3, 3, Float64, 9}(I), R::Float64 = 1.0)
    nside >= 1 || throw(ArgumentError("nside must be >= 1, got $nside"))
    ordering in (:ring, :nested) ||
        throw(ArgumentError("ordering must be :ring or :nested, got $ordering"))
    R > 0 || throw(ArgumentError("R must be positive, got $R"))
    rotation = convert(SMatrix{3, 3, Float64, 9}, rotation)

    M = Sphere(2)
    n_cells = 12 * nside * nside

    n_rings = 4 * nside - 1

    # --- Step 1: Generate cell centers on rings ---
    # ring_centers[ring][i] = center of cell i in ring (1-indexed)
    ring_centers = Vector{Vector{SVector{3, Float64}}}(undef, n_rings)
    ring_counts = Vector{Int}(undef, n_rings)

    for ring in 1:n_rings
        if ring <= nside
            # North polar cap
            n_in_ring = 4 * ring
            cos_theta = 1.0 - ring^2 / (3.0 * nside^2)
            theta = acos(cos_theta)
        elseif ring <= 3 * nside
            # Equatorial belt
            n_in_ring = 4 * nside
            cos_theta = (4.0 * nside - 2.0 * ring) / (3.0 * nside)
            theta = acos(cos_theta)
        else
            # South polar cap
            j = n_rings - ring + 1
            n_in_ring = 4 * j
            cos_theta = -(1.0 - j^2 / (3.0 * nside^2))
            theta = acos(cos_theta)
        end

        ring_counts[ring] = n_in_ring
        centers = SVector{3, Float64}[]

        for i in 1:n_in_ring
            phi = 2π * (i - 0.5) / n_in_ring
            x = sin(theta) * cos(phi)
            y = sin(theta) * sin(phi)
            z = cos(theta)
            p = rotation * SVector(x, y, z)
            push!(centers, SVector{3, Float64}(p))
        end
        ring_centers[ring] = centers
    end

    # --- Step 2: Compute corner vertices for each cell ---
    # For each cell, the 4 corners are at the intersections of ring boundaries
    # and sector boundaries. Each corner is the normalized average of the
    # adjacent cell centers that surround that corner.
    # For boundary corners at poles, use the pole point + adjacent centers.

    # Helper to get a ring center with periodic wrap
    function get_center(ring, i)
        n = ring_counts[ring]
        idx = mod(i - 1, n) + 1
        return ring_centers[ring][idx]
    end

    north_pole = rotation * SVector(0.0, 0.0, 1.0)
    south_pole = rotation * SVector(0.0, 0.0, -1.0)

    # Deduplicate corners using a Dict with rounded coordinates as keys
    corner_dict = Dict{Tuple{Int, Int, Int}, Int}()
    nodes = SVector{3, Float64}[]
    _cell_nodes = Vector{NTuple{4, Int}}(undef, n_cells)

    function add_corner(v::SVector{3, Float64})
        # Round to 12 decimal places for deduplication (~1e-12 m precision at R=1)
        key = (round(Int, v[1] * 1e12), round(Int, v[2] * 1e12), round(Int, v[3] * 1e12))
        if haskey(corner_dict, key)
            return corner_dict[key]
        end
        idx = length(nodes) + 1
        push!(nodes, v)
        corner_dict[key] = idx
        return idx
    end

    function normalized_mean(vectors)
        s = sum(vectors)
        n = norm(s)
        if n < 1e-15
            error("normalized_mean: zero vector")
        end
        return R * (s / n)
    end

    cell_id = 1
    for ring in 1:n_rings
        n_in_ring = ring_counts[ring]
        for i in 1:n_in_ring
            c = get_center(ring, i)
            w = get_center(ring, i - 1)
            e = get_center(ring, i + 1)

            # Determine north and south neighbors
            if ring > 1
                n_ring = ring - 1
                n_count = ring_counts[n_ring]
                # Find the cell in the north ring that is closest in longitude
                # The north ring has fewer or equal cells
                # For HEALPix ring ordering, cell i in ring maps to approximately
                # the same longitudinal sector in the north ring
                ratio = n_count / n_in_ring
                i_north = clamp(round(Int, (i - 0.5) * ratio + 0.5), 1, n_count)
                n = get_center(n_ring, i_north)
                nw = get_center(n_ring, i_north - 1)
                ne = get_center(n_ring, i_north + 1)
            else
                # North pole
                n = nothing
            end

            if ring < n_rings
                s_ring = ring + 1
                s_count = ring_counts[s_ring]
                ratio = s_count / n_in_ring
                i_south = clamp(round(Int, (i - 0.5) * ratio + 0.5), 1, s_count)
                s = get_center(s_ring, i_south)
                sw_s = get_center(s_ring, i_south - 1)
                se_s = get_center(s_ring, i_south + 1)
            else
                # South pole
                s = nothing
            end

            # Compute 4 corners (SW, SE, NE, NW in local cell coordinates)
            # For ring 1 (north polar cap), the "north" side is the pole.
            # For ring n_rings (south polar cap), the "south" side is the pole.

            # SW corner: on the southern boundary, western side
            if ring == n_rings
                # South polar cap: SW corner is the south pole
                sw = south_pole
            else
                sw = normalized_mean([c, w, s, sw_s])
            end

            # SE corner: on the southern boundary, eastern side
            if ring == n_rings
                # South polar cap: SE corner is the south pole
                se = south_pole
            else
                se = normalized_mean([c, e, s, se_s])
            end

            # NW corner: on the northern boundary, western side
            if ring == 1
                # North polar cap: NW corner is the north pole
                nw_corner = north_pole
            else
                nw_corner = normalized_mean([c, w, n, nw])
            end

            # NE corner: on the northern boundary, eastern side
            if ring == 1
                # North polar cap: NE corner is the north pole
                ne_corner = north_pole
            else
                ne_corner = normalized_mean([c, e, n, ne])
            end

            sw_id = add_corner(sw)
            se_id = add_corner(se)
            ne_id = add_corner(ne_corner)
            nw_id = add_corner(nw_corner)

            _cell_nodes[cell_id] = (sw_id, se_id, ne_id, nw_id)
            cell_id += 1
        end
    end

    # --- Step 3: Compute cell volumes and centroids from corner vertices ---
    cell_volumes = Vector{Float64}(undef, n_cells)
    cell_centroids = Vector{SVector{3, Float64}}(undef, n_cells)

    for cid in 1:n_cells
        cn = _cell_nodes[cid]
        A = nodes[cn[1]]   # SW
        B = nodes[cn[2]]   # SE
        C = nodes[cn[3]]   # NE
        D = nodes[cn[4]]   # NW

        area = spherical_triangle_area(R, A, B, C) +
               spherical_triangle_area(R, A, C, D)
        cell_volumes[cid] = area

        verts = [A, B, C, D]
        c = Manifolds.mean(M, verts)
        cell_centroids[cid] = SVector{3, Float64}(c)
    end

    # Scale volumes to enforce exact area conservation (compensates for
    # non-conforming gaps/overlaps in the approximate corner reconstruction).
    total_area = sum(cell_volumes)
    if total_area > 0
        scale = 4π * R^2 / total_area
        for cid in 1:n_cells
            cell_volumes[cid] *= scale
        end
    end

    # --- Derive edges from cell-node connectivity ---
    edge_map = Dict{Tuple{Int, Int}, Int}()
    _cell_edges = NTuple{4, Int}[]

    for cell_id in 1:n_cells
        cn = _cell_nodes[cell_id]
        cell_edge_ids = Int[]
        for (a, b) in ((cn[1], cn[2]), (cn[2], cn[3]),
                        (cn[3], cn[4]), (cn[4], cn[1]))
            key = a < b ? (a, b) : (b, a)
            edge_id = get!(edge_map, key) do
                length(edge_map) + 1
            end
            push!(cell_edge_ids, edge_id)
        end
        push!(_cell_edges, tuple(cell_edge_ids...))
    end

    n_edges = length(edge_map)
    _edge_nodes = Vector{NTuple{2, Int}}(undef, n_edges)
    for ((n1, n2), edge_id) in edge_map
        _edge_nodes[edge_id] = (n1, n2)
    end

    # --- Derive cell neighbors from edge sharing ---
    edge_cells = [Int[] for _ in 1:n_edges]
    for cell_id in 1:n_cells
        for e in _cell_edges[cell_id]
            push!(edge_cells[e], cell_id)
        end
    end

    _cell_cells = Vector{NTuple{4, Int}}(undef, n_cells)
    for cell_id in 1:n_cells
        ce = _cell_edges[cell_id]
        neighbors = Int[]
        for e in ce
            adj = edge_cells[e]
            if length(adj) == 1
                # Boundary edge or self-loop — no neighbor
                push!(neighbors, 0)
            else
                other = first(c for c in adj if c != cell_id)
                push!(neighbors, other)
            end
        end
        _cell_cells[cell_id] = tuple(neighbors...)
    end

    # --- Apply nested ordering permutation if requested ---
    if ordering == :nested
        # Flatten ring centers for permutation computation
        ring_centers_flat = SVector{3, Float64}[]
        for ring in 1:n_rings
            append!(ring_centers_flat, ring_centers[ring])
        end
        perm = _ring_to_nested_permutation(nside, ring_centers_flat)
        # perm[nested_id] = ring_id

        # Build inverse permutation: inv_perm[ring_id] = nested_id
        inv_perm = Vector{Int}(undef, n_cells)
        for nested_id in 1:n_cells
            inv_perm[perm[nested_id]] = nested_id
        end

        # Apply permutation to cell-related arrays
        cell_volumes = [cell_volumes[perm[i]] for i in 1:n_cells]
        cell_centroids = [cell_centroids[perm[i]] for i in 1:n_cells]
        _cell_nodes = [_cell_nodes[perm[i]] for i in 1:n_cells]
        _cell_edges = [_cell_edges[perm[i]] for i in 1:n_cells]

        # For cell neighbors, permute the neighbor IDs too
        _cell_cells_perm = Vector{NTuple{4, Int}}(undef, n_cells)
        for i in 1:n_cells
            old_neighbors = _cell_cells[perm[i]]
            new_neighbors = ntuple(j -> old_neighbors[j] == 0 ? 0 : inv_perm[old_neighbors[j]], 4)
            _cell_cells_perm[i] = new_neighbors
        end
        _cell_cells = _cell_cells_perm
    end

    return HEALPixGrid{typeof(M)}(
        M, nside, R, ordering, nodes,
        cell_volumes, cell_centroids, _cell_nodes,
        _cell_edges, _cell_cells, _edge_nodes,
        Ref{Union{Nothing, AbstractManifoldMesh{typeof(M)}}}(nothing),
        rotation)
end

# -- Trait Implementations --

TopologyStyle(::Type{<:HEALPixGrid}) = IsSemiGrid()
CellTypeStyle(::Type{<:HEALPixGrid}) = IsUniform{4}()
PatchStyle(::Type{<:HEALPixGrid}) = NoPatch()

has_dual(g::HEALPixGrid) = g._dual[] !== nothing

# -- Global Information --

manifold(g::HEALPixGrid) = g.manifold
num_cells(g::HEALPixGrid) = 12 * g.nside * g.nside
num_nodes(g::HEALPixGrid) = length(g.nodes)
num_edges(g::HEALPixGrid) = length(g._edge_nodes)

# -- Geometry --

function node_coordinates(g::HEALPixGrid, node_id::Int)
    _check_node_id(g, node_id)
    return g.nodes[node_id]
end

function cell_volume(g::HEALPixGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g.cell_volumes[cell_id]
end

function cell_centroid(g::HEALPixGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g.cell_centroids[cell_id]
end

# -- Topology --

function cell_nodes(g::HEALPixGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g._cell_nodes[cell_id]
end

function cell_cells(g::HEALPixGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g._cell_cells[cell_id]
end

function node_cells(g::HEALPixGrid, node_id::Int)
    _check_node_id(g, node_id)
    cells = Int[]
    for cell_id in 1:num_cells(g)
        if node_id in cell_nodes(g, cell_id)
            push!(cells, cell_id)
        end
    end
    return cells
end

function cell_edges(g::HEALPixGrid, cell_id::Int)
    _check_cell_id(g, cell_id)
    return g._cell_edges[cell_id]
end

# -- Edge Geometry --

@inline function _check_edge_id(g::HEALPixGrid, edge_id::Int)
    @boundscheck 1 <= edge_id <= num_edges(g) ||
                 throw(BoundsError("edge_id $edge_id out of range [1, $(num_edges(g))]"))
    nothing
end

function _edge_endpoints(g::HEALPixGrid, edge_id::Int)
    n1, n2 = g._edge_nodes[edge_id]
    return (node_coordinates(g, n1), node_coordinates(g, n2))
end

function edge_length(g::HEALPixGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    n1, n2 = _edge_endpoints(g, edge_id)
    return Manifolds.distance(g.manifold, n1, n2)
end

function edge_midpoint(g::HEALPixGrid, edge_id::Int)
    _check_edge_id(g, edge_id)
    n1, n2 = _edge_endpoints(g, edge_id)
    if Manifolds.distance(g.manifold, n1, n2) < 1e-14
        return SVector{3, Float64}(n1)
    end
    return SVector{3, Float64}(Manifolds.mid_point(g.manifold, n1, n2))
end

function edge_outward_normal(g::HEALPixGrid, edge_id::Int, cell_id::Int)
    _check_edge_id(g, edge_id)
    _check_cell_id(g, cell_id)
    n1, n2 = _edge_endpoints(g, edge_id)

    if Manifolds.distance(g.manifold, n1, n2) < 1e-14
        midpoint = n1
        return (base_point = SVector{3, Float64}(midpoint),
            normal = zero(SVector{3, Float64}))
    end

    midpoint = Manifolds.mid_point(g.manifold, n1, n2)
    gc_normal = cross(SVector(n1), SVector(n2))
    tangent = normalize(cross(gc_normal, SVector(midpoint)))
    cell_c = cell_centroid(g, cell_id)
    cell_side = sign(dot(gc_normal, SVector(cell_c)))
    outward = cell_side * cross(tangent, SVector(midpoint))
    outward = Manifolds.project(g.manifold, midpoint, outward)

    return (base_point = SVector{3, Float64}(midpoint),
        normal = SVector{3, Float64}(outward))
end

# -- Boundary --

boundary_nodes(g::HEALPixGrid, marker) = Int[]
boundary_edges(g::HEALPixGrid, marker) = Int[]
