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
              view::Symbol = Symbol("3d"),
              projection::Symbol = :equirectangular,
              n_arc_points::Int = 20,
              figsize::Tuple{Int,Int} = (800, 600),
              kwargs...) -> Figure

Plot mesh wireframe. Returns a Makie `Figure`.
"""
function plot_mesh(g::AbstractManifoldMesh;
                   show_nodes::Bool = false,
                   show_edges::Bool = true,
                   show_cell_ids::Bool = false,
                   show_node_ids::Bool = false,
                   view::Symbol = Symbol("3d"),
                   projection::Symbol = :equirectangular,
                   n_arc_points::Int = 20,
                   figsize::Tuple{Int,Int} = (800, 600),
                   kwargs...)
    view == Symbol("3d") ? _plot_mesh_3d(g; show_nodes, show_edges, show_cell_ids,
                                  show_node_ids, n_arc_points, figsize, kwargs...) :
    view == Symbol("2d") ? _plot_mesh_2d(g; show_nodes, show_edges, show_cell_ids,
                                  show_node_ids, n_arc_points, projection,
                                  figsize, kwargs...) :
    error("Unsupported view: $view. Use Symbol(\"3d\") or Symbol(\"2d\").")
end

"""
    plot_mesh_filled(g::AbstractManifoldMesh;
                     show_edges::Bool = true,
                     color_by = nothing,
                     kwargs...) -> Figure

Plot filled cell polygons with optional wireframe overlay. Returns a Makie `Figure`.

# Keywords
- `color_by = nothing`: callable taking cell index (`Int`) and returning a Makie color,
  or `nothing` for default `lightblue`.
"""
function plot_mesh_filled(g::AbstractManifoldMesh;
                          show_edges::Bool = true,
                          color_by = nothing,
                          view::Symbol = Symbol("3d"),
                          n_arc_points::Int = 20,
                          figsize::Tuple{Int,Int} = (800, 600),
                          kwargs...)
    view == Symbol("3d") ? _plot_filled_3d(g; show_edges, color_by, n_arc_points, figsize, kwargs...) :
    view == Symbol("2d") ? _plot_filled_2d(g; show_edges, color_by, n_arc_points, figsize, kwargs...) :
    error("Unsupported view: $view. Use Symbol(\"3d\") or Symbol(\"2d\").")
end

# -- 3D wireframe --

function _plot_mesh_3d(g; show_nodes, show_edges, show_cell_ids,
                        show_node_ids, n_arc_points, figsize, kwargs...)
    M = _require_makie()

    fig = M.Figure(; size=figsize)
    ax = M.LScene(fig[1, 1]; show_axis=false)

    _add_sphere_background!(ax, M; R=_get_radius(g))

    if show_edges
        segs = edge_segments(g; n_arc_points)
        for seg in segs
            if length(seg) > 1
                M.lines!(ax, seg; color=:steelblue, linewidth=0.5, kwargs...)
            end
        end
    end

    if show_nodes
        pts = node_points(g)
        M.meshscatter!(ax, pts; color=:orangered, markersize=5, kwargs...)
    end

    if show_cell_ids
        for cid in 1:num_cells(g)
            c = cell_centroid(g, cid)
            M.text!(ax, [M.Point3f(Float32.(c))],
                    text=["$cid"], fontsize=6, color=:black, kwargs...)
        end
    end

    if show_node_ids
        for nid in 1:num_nodes(g)
            c = node_coordinates(g, nid)
            M.text!(ax, [M.Point3f(Float32.(c))],
                    text=["$nid"], fontsize=5, color=:gray, kwargs...)
        end
    end

    return fig
end

# -- 3D filled --

function _plot_filled_3d(g; show_edges, color_by, n_arc_points, figsize, kwargs...)
    M = _require_makie()

    fig = M.Figure(; size=figsize)
    ax = M.LScene(fig[1, 1]; show_axis=false)

    _add_sphere_background!(ax, M; R=_get_radius(g))

    polys = cell_polygons(g; n_arc_points)
    for (i, poly) in enumerate(polys)
        if length(poly) >= 3
            center = poly[1]
            color = color_by === nothing ? :lightblue : color_by(i)
            for k in 2:length(poly)-1
                M.mesh!(ax, [center, poly[k], poly[k+1]];
                        color=color, transparency=true, kwargs...)
            end
        end
    end

    if show_edges
        segs = edge_segments(g; n_arc_points)
        for seg in segs
            if length(seg) > 1
                M.lines!(ax, seg; color=:steelblue, linewidth=0.5, kwargs...)
            end
        end
    end

    return fig
end

# -- 2D wireframe (equirectangular) --

function _plot_mesh_2d(g; show_nodes, show_edges, show_cell_ids,
                        show_node_ids, n_arc_points, projection,
                        figsize, kwargs...)
    projection == :equirectangular || error("Only :equirectangular projection is supported")
    M = _require_makie()

    fig = M.Figure(; size=figsize)
    ax = M.Axis(fig[1, 1];
                xlabel="Longitude (deg)", ylabel="Latitude (deg)",
                title="Mesh (equirectangular)")

    if show_edges
        segs = edge_segments(g; n_arc_points)
        for seg in segs
            if length(seg) > 1
                lons = [_lon_from_xyz(p) for p in seg]
                lats = [_lat_from_xyz(p) for p in seg]
                M.lines!(ax, lons, lats; color=:steelblue, linewidth=0.5, kwargs...)
            end
        end
    end

    if show_nodes
        pts = node_points(g)
        xs = [_lon_from_xyz(p) for p in pts]
        ys = [_lat_from_xyz(p) for p in pts]
        M.scatter!(ax, xs, ys; color=:orangered, markersize=5, kwargs...)
    end

    if show_cell_ids
        for cid in 1:num_cells(g)
            c = cell_centroid(g, cid)
            lon = _lon_from_xyz(Point3f(Float32.(c)))
            lat = _lat_from_xyz(Point3f(Float32.(c)))
            M.text!(ax, [lon], [lat]; text=["$cid"], fontsize=6, kwargs...)
        end
    end

    if show_node_ids
        for nid in 1:num_nodes(g)
            c = node_coordinates(g, nid)
            lon = _lon_from_xyz(Point3f(Float32.(c)))
            lat = _lat_from_xyz(Point3f(Float32.(c)))
            M.text!(ax, [lon], [lat]; text=["$nid"], fontsize=5, color=:gray, kwargs...)
        end
    end

    return fig
end

# -- 2D filled --

function _plot_filled_2d(g; show_edges, color_by, n_arc_points, figsize, kwargs...)
    M = _require_makie()

    fig = M.Figure(; size=figsize)
    ax = M.Axis(fig[1, 1];
                xlabel="Longitude (deg)", ylabel="Latitude (deg)",
                title="Mesh (equirectangular)")

    polys = cell_polygons(g; n_arc_points)
    for poly in polys
        if length(poly) >= 3
            lons = [_lon_from_xyz(p) for p in poly]
            lats = [_lat_from_xyz(p) for p in poly]
            push!(lons, lons[1])
            push!(lats, lats[1])
            M.lines!(ax, lons, lats; color=:lightblue, linewidth=1, kwargs...)
            M.poly!(ax, lons, lats; color=:lightblue, strokewidth=0, kwargs...)
        end
    end

    if show_edges
        segs = edge_segments(g; n_arc_points)
        for seg in segs
            if length(seg) > 1
                lons = [_lon_from_xyz(p) for p in seg]
                lats = [_lat_from_xyz(p) for p in seg]
                M.lines!(ax, lons, lats; color=:steelblue, linewidth=0.5, kwargs...)
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
    return sqrt(sum(c.^2))
end

function _add_sphere_background!(ax, M; R=1.0)
    n = 64
    theta = LinRange(0, pi, n)
    phi = LinRange(-pi, pi, 2n)
    xe = [R * cos(phiv) * sin(thetav) for thetav in theta, phiv in phi]
    ye = [R * sin(phiv) * sin(thetav) for thetav in theta, phiv in phi]
    ze = [R * cos(thetav) for thetav in theta, phiv in phi]
    M.surface!(ax, xe, ye, ze;
              color=(:lightgray, 0.15),
              transparency=true,
              shading=M.NoShading)
end

_lon_from_xyz(p) = rad2deg(atan(p[2], p[1]))
_lat_from_xyz(p) = rad2deg(asin(clamp(p[3] / sqrt(sum(p.^2)), -1, 1)))

end  # module VisualizationPlotting
