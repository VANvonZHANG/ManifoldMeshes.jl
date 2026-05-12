@testset "TopologyStyle trait" begin
    @test TopologyStyle(IsGrid) === IsGrid()
    @test TopologyStyle(IsMesh) === IsMesh()
    @test isa(TopologyStyle(IsGrid), TopologyStyle)
end

@testset "AbstractLocation types" begin
    @test NodeLoc <: AbstractLocation
    @test CellLoc <: AbstractLocation
    @test EdgeLoc <: AbstractLocation
end

@testset "TopologyStyle for LatLonGrid" begin
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 90.0], lon_edges = [0.0, 360.0])
    @test TopologyStyle(g) === IsGrid()
    @test TopologyStyle(LatLonGrid) === IsGrid()
end

@testset "MixedCellTopology" begin
    m = ManifoldMeshes.MixedCellTopology((1, 2, 3, 0, 0, 0), 3)
    @test length(m) == 3
    @test size(m) == (3,)
    @test m[1] == 1
    @test m[2] == 2
    @test m[3] == 3
    @test collect(m) == [1, 2, 3]
    @test_throws BoundsError m[4]
end
