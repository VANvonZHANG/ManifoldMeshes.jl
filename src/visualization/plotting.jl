module VisualizationPlotting

using ..ManifoldMeshes
using ..ManifoldMeshes: node_points, edge_segments, cell_polygons
using GeometryBasics: Point3f

"""
    plot_mesh(g::AbstractManifoldMesh;
              show_nodes::Bool = false,
              show_edges::Bool = true,
              show_cell_ids::Bool = false,
              show_node_ids::Bool = false,
              n_arc_points::Int = 50,
              figsize::Tuple{Int,Int} = (800, 600),
              kwargs...) -> Figure

Plot 3D mesh wireframe on a semi-transparent sphere background.
Returns a Makie `Figure`.
"""
function plot_mesh(g::AbstractManifoldMesh;
        show_nodes::Bool = false,
        show_edges::Bool = true,
        show_cell_ids::Bool = false,
        show_node_ids::Bool = false,
        n_arc_points::Int = 50,
        figsize::Tuple{Int, Int} = (800, 600),
        kwargs...)
    M = _require_makie()

    fig = M.Figure(; size = figsize)
    ax = M.LScene(fig[1, 1]; show_axis = false)

    _add_sphere_background!(ax, M; R = _get_radius(g))

    if show_edges
        segs = edge_segments(g; n_arc_points)
        for seg in segs
            if length(seg) > 1
                M.lines!(ax, seg; color = :steelblue, linewidth = 0.5, kwargs...)
            end
        end
    end

    if show_nodes
        pts = node_points(g)
        M.meshscatter!(ax, pts; color = :orangered, markersize = 0.03, kwargs...)
    end

    if show_cell_ids
        for cid in 1:num_cells(g)
            c = cell_centroid(g, cid)
            M.text!(ax, [Point3f(Float32.(c))],
                text = ["$cid"], fontsize = 6, color = :black, kwargs...)
        end
    end

    if show_node_ids
        for nid in 1:num_nodes(g)
            c = node_coordinates(g, nid)
            M.text!(ax, [Point3f(Float32.(c))],
                text = ["$nid"], fontsize = 5, color = :gray, kwargs...)
        end
    end

    return fig
end

"""
    plot_mesh_filled(g::AbstractManifoldMesh;
                     show_edges::Bool = true,
                     color_by = nothing,
                     n_arc_points::Int = 50,
                     figsize::Tuple{Int,Int} = (800, 600),
                     kwargs...) -> Figure

Plot filled 3D cell polygons on a semi-transparent sphere background.
Returns a Makie `Figure`.

# Keywords
- `color_by = nothing`: callable taking cell index (`Int`) and returning a Makie color,
  or `nothing` for default `lightblue`.
"""
function plot_mesh_filled(g::AbstractManifoldMesh;
        show_edges::Bool = true,
        color_by = nothing,
        n_arc_points::Int = 50,
        figsize::Tuple{Int, Int} = (800, 600),
        kwargs...)
    M = _require_makie()

    fig = M.Figure(; size = figsize)
    ax = M.LScene(fig[1, 1]; show_axis = false)

    _add_sphere_background!(ax, M; R = _get_radius(g))

    polys = cell_polygons(g; n_arc_points)
    for (i, poly) in enumerate(polys)
        if length(poly) >= 3
            center = poly[1]
            color = color_by === nothing ? :lightblue : color_by(i)
            for k in 2:(length(poly) - 1)
                M.mesh!(ax, [center, poly[k], poly[k + 1]];
                    color = color, transparency = true, shading = M.NoShading, kwargs...)
            end
        end
    end

    if show_edges
        segs = edge_segments(g; n_arc_points)
        for seg in segs
            if length(seg) > 1
                M.lines!(ax, seg; color = :steelblue, linewidth = 0.5, kwargs...)
            end
        end
    end

    return fig
end

# -- Helpers --

function _require_makie()
    if !isdefined(Main, :Makie)
        error("Makie is required for visualization. Run: `using CairoMakie` or `using GLMakie` before calling plot_mesh.")
    end
    return Main.Makie
end

function _get_radius(g::AbstractManifoldMesh)
    c = node_coordinates(g, 1)
    return sqrt(sum(c .^ 2))
end

function _add_sphere_background!(ax, M; R = 1.0)
    n = 64
    theta = LinRange(0, pi, n)
    phi = LinRange(-pi, pi, 2n)
    xe = [R * cos(phiv) * sin(thetav) for thetav in theta, phiv in phi]
    ye = [R * sin(phiv) * sin(thetav) for thetav in theta, phiv in phi]
    ze = [R * cos(thetav) for thetav in theta, phiv in phi]
    colors = fill(M.RGBA(0.9, 0.9, 0.9, 0.15), size(xe))
    M.surface!(ax, xe, ye, ze;
        color = colors,
        transparency = true,
        shading = M.NoShading)
end

end  # module VisualizationPlotting
