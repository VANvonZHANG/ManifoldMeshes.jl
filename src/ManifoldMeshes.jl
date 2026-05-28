module ManifoldMeshes

using Manifolds
using ManifoldsBase
using StaticArrays
using LinearAlgebra
using GeometryBasics: Point3f
using CairoMakie

include("traits.jl")
include("interface.jl")
include("dual.jl")
include("csr.jl")

# Export types and functions
export TopologyStyle, IsGrid, IsSemiGrid, IsMesh, AbstractLocation, NodeLoc, CellLoc,
       EdgeLoc
export CellTypeStyle, IsUniform, IsMixed
export PatchStyle, NoPatch, MultiPatch
export AbstractManifoldMesh, AbstractDualMesh, LatLonGrid, CubedSphereGrid,
       ReducedGaussianGrid, HEALPixGrid, MixedCellTopology
export manifold, num_cells, num_nodes, num_edges
export node_coordinates, cell_volume, cell_centroid
export cell_nodes, cell_cells, node_cells, cell_edges
export edge_length, edge_midpoint, edge_outward_normal
export boundary_nodes, boundary_edges
export dual, has_dual, cell_face, cell_local_2d
export slerp, node_points, edge_segments, cell_polygons, cell_triangles

include("sphere/utils.jl")
include("sphere/latlon.jl")
include("sphere/cubed_sphere.jl")
include("sphere/reduced_gaussian.jl")
include("sphere/healpix.jl")
include("visualization/mesh_data.jl")
include("visualization/plotting.jl")

for sym in [:plot_mesh, :plot_mesh_filled]
    @eval ManifoldMeshes const $sym = VisualizationPlotting.$sym
end

export plot_mesh, plot_mesh_filled

end # module
