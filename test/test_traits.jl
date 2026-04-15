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