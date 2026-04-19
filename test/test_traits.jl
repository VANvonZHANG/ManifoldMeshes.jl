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
