# examples/grids_visualization.jl
# Standalone script: visualize all 4 S² grid types at coarse and fine resolutions.
#
# Run with:
#   julia --project=. examples/grids_visualization.jl
#
# Output:
#   examples/grids_visualization.png

using ManifoldMeshes
using CairoMakie
using GeometryBasics: Point3f
using LinearAlgebra: norm

# ── Tunable constants ──────────────────────────────────────────────────────────

const OUTPUT_PATH = joinpath(@__DIR__, "grids_visualization.png")
const FIG_SIZE = (900, 1100)    # (width, height) in pixels
const N_ARC = 30             # edge discretization points

# ── Grid definitions (coarse / fine) ─────────────────────────────────────────

# LatLon: structured lat-lon grid
#   Coarse: 6 lat bands × 8 lon bands  = 48 cells
#   Fine:   24 lat bands × 32 lon bands = 768 cells
latlon_coarse = LatLonGrid(
    lat_edges = collect(range(-90.0, 90.0, length = 7)),
    lon_edges = collect(range(0.0, 360.0, length = 9))
)
latlon_fine = LatLonGrid(
    lat_edges = collect(range(-90.0, 90.0, length = 25)),
    lon_edges = collect(range(0.0, 360.0, length = 33))
)

# CubedSphere: 6-face gnomonic projection
#   Coarse: n=4  → 6 × 4²  = 96 cells
#   Fine:   n=16 → 6 × 16² = 1536 cells
cubed_coarse = CubedSphereGrid(n = 4)
cubed_fine = CubedSphereGrid(n = 16)

# ReducedGaussian: octahedral Gaussian lat bands
#   Coarse: nlat=8  → ~96 cells
#   Fine:   nlat=32 → ~1536 cells
gauss_coarse = ReducedGaussianGrid(nlat = 8)
gauss_fine = ReducedGaussianGrid(nlat = 32)

# HEALPix: hierarchical equal-area pixels
#   Coarse: nside=2 → 12 × 2²  = 48 cells
#   Fine:   nside=8 → 12 × 8²  = 768 cells
healpix_coarse = HEALPixGrid(nside = 2)
healpix_fine = HEALPixGrid(nside = 8)

# ── Bundle for iteration ──────────────────────────────────────────────────────

# (name, coarse_grid, fine_grid, color)
grids = [
    ("LatLon", latlon_coarse, latlon_fine, :steelblue3),
    ("CubedSphere", cubed_coarse, cubed_fine, :seagreen),
    ("ReducedGaussian", gauss_coarse, gauss_fine, :mediumpurple3),
    ("HEALPix", healpix_coarse, healpix_fine, :chocolate3)
]

# ── Helper: draw mesh wireframe on an existing LScene ─────────────────────────

function draw_mesh!(ax, g; edge_color = :steelblue, show_edges::Bool = true)
    # Semi-transparent sphere background
    n = 64
    theta = LinRange(0, pi, n)
    phi = LinRange(-pi, pi, 2n)
    R = norm(node_coordinates(g, 1))
    xe = [R * cos(phiv) * sin(thetav) for thetav in theta, phiv in phi]
    ye = [R * sin(phiv) * sin(thetav) for thetav in theta, phiv in phi]
    ze = [R * cos(thetav) for thetav in theta, phiv in phi]
    colors = fill(RGBAf(0.85, 0.85, 0.9, 0.25), size(xe))
    surface!(ax, xe, ye, ze;
        color = colors, transparency = true, shading = NoShading)

    if show_edges
        segs = edge_segments(g; n_arc_points = N_ARC)
        pts = Point3f[]
        nan_pt = Point3f(NaN32, NaN32, NaN32)
        for seg in segs
            if length(seg) >= 2
                append!(pts, seg)
                push!(pts, nan_pt)
            end
        end
        isempty(pts) || lines!(ax, pts; color = edge_color, linewidth = 0.8)
    end

    return nothing
end

# ── Create figure and plot ────────────────────────────────────────────────────

fig = Figure(size = FIG_SIZE)

for (row, (name, coarse, fine, color)) in enumerate(grids)
    # Coarse resolution (left column)
    Label(fig[row, 1, Top()], "$name (coarse, $(num_cells(coarse)) cells)";
        fontsize = 12, font = :bold, padding = (0, 0, 5, 0))
    ax1 = LScene(fig[row, 1]; show_axis = false)
    draw_mesh!(ax1, coarse; edge_color = color)

    # Fine resolution (right column)
    Label(fig[row, 2, Top()], "$name (fine, $(num_cells(fine)) cells)";
        fontsize = 12, font = :bold, padding = (0, 0, 5, 0))
    ax2 = LScene(fig[row, 2]; show_axis = false)
    draw_mesh!(ax2, fine; edge_color = color)
end

# Tighten layout gaps
rowgap!(fig.layout, 15)
colgap!(fig.layout, 25)

# ── Save ──────────────────────────────────────────────────────────────────────

save(OUTPUT_PATH, fig)
println("Saved visualization to: $OUTPUT_PATH")
