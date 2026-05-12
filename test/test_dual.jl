using ManifoldMeshes
using Test

@testset "Dual framework stubs" begin
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 90.0], lon_edges = [0.0, 360.0])
    @test has_dual(g) == false
    # LatLonGrid does not yet support dual construction
    @test_throws ErrorException dual(g)
end
