using ManifoldMeshes
using Manifolds
using StaticArrays
using Test

@testset "ManifoldMeshes.jl" begin
    include("test_traits.jl")
    include("test_latlon_construction.jl")
    include("test_latlon_geometry.jl")
    include("test_latlon_connectivity.jl")
    include("test_latlon_normals.jl")
end
