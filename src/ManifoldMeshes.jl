module ManifoldMeshes

using Manifolds
using ManifoldsBase
using StaticArrays
using LinearAlgebra

include("traits.jl")
include("interface.jl")

# Export types and functions
export TopologyStyle, IsGrid, IsMesh, AbstractLocation, NodeLoc, CellLoc, EdgeLoc
export AbstractManifoldMesh
include("sphere/latlon.jl")

end # module