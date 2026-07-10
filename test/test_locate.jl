using ManifoldMeshes
using StaticArrays
using Test

@testset "locate helpers" begin
    # _bilinear_weights: corners are unit vectors
    @test ManifoldMeshes._bilinear_weights(0.0, 0.0) == (1.0, 0.0, 0.0, 0.0)   # SW
    @test ManifoldMeshes._bilinear_weights(1.0, 0.0) == (0.0, 1.0, 0.0, 0.0)   # SE
    @test ManifoldMeshes._bilinear_weights(1.0, 1.0) == (0.0, 0.0, 1.0, 0.0)   # NE
    @test ManifoldMeshes._bilinear_weights(0.0, 1.0) == (0.0, 0.0, 0.0, 1.0)   # NW
    # centroid: all 0.25
    @test ManifoldMeshes._bilinear_weights(0.5, 0.5) == (0.25, 0.25, 0.25, 0.25)
    # partition of unity
    for (s, t) in [(0.3, 0.7), (0.9, 0.1), (0.123, 0.456)]
        @test sum(ManifoldMeshes._bilinear_weights(s, t)) ≈ 1.0
    end
end

@testset "coordinate round-trip" begin
    for (lat, lon) in
        [(0.0, 0.0), (45.0, 90.0), (-30.0, 200.0), (89.9, -175.0), (0.0, 359.9)]
        p = ManifoldMeshes._latlon_to_cartesian(lat, lon, 1.0)
        lat2, lon2 = ManifoldMeshes._cartesian_to_latlon(p)
        @test lat2 ≈ lat atol = 1e-9
        @test min(mod(lon2 - lon, 360.0), 360.0 - mod(lon2 - lon, 360.0)) ≈ 0 atol = 1e-9
    end
    # pole: longitude is arbitrary but latitude exact
    p = ManifoldMeshes._latlon_to_cartesian(90.0, 42.0)
    lat2, lon2 = ManifoldMeshes._cartesian_to_latlon(p)
    @test lat2 ≈ 90.0 atol = 1e-9
end
