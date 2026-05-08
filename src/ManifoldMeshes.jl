module ManifoldMeshes

using Manifolds
using ManifoldsBase
using StaticArrays
using LinearAlgebra
using GeometryBasics: Point3f

include("traits.jl")
include("interface.jl")

# Export types and functions
export TopologyStyle, IsGrid, IsMesh, AbstractLocation, NodeLoc, CellLoc, EdgeLoc
export AbstractManifoldMesh, LatLonGrid
export manifold, num_cells, num_nodes, num_edges
export node_coordinates, cell_volume, cell_centroid
export cell_nodes, cell_cells, node_cells, cell_edges
export edge_length, edge_midpoint, edge_outward_normal
export boundary_nodes, boundary_edges
export slerp, node_points, edge_segments, cell_polygons

include("sphere/latlon.jl")
include("visualization/mesh_data.jl")

end # module
