using ManifoldMeshes
using Manifolds
using StaticArrays
using Test
using Aqua

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(ManifoldMeshes)
end

@testset "ManifoldMeshes.jl" begin
    include("test_traits.jl")
    include("test_latlon_construction.jl")
    include("test_latlon_geometry.jl")
    include("test_latlon_connectivity.jl")
    include("test_latlon_normals.jl")
    include("test_latlon_edge_cases.jl")
end
