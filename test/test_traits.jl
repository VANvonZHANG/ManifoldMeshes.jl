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
    # basic usage
    m = ManifoldMeshes.MixedCellTopology((1, 2, 3, 0, 0, 0), 3)
    @test length(m) == 3
    @test size(m) == (3,)
    @test m[1] == 1
    @test m[2] == 2
    @test m[3] == 3
    @test collect(m) == [1, 2, 3]
    @test_throws BoundsError m[4]

    # edge case: empty topology (len = 0)
    m0 = ManifoldMeshes.MixedCellTopology((0, 0, 0), 0)
    @test length(m0) == 0
    @test size(m0) == (0,)
    @test collect(m0) == Int[]
    @test_throws BoundsError m0[1]

    # edge case: full topology (len = MAX_K)
    m_full = ManifoldMeshes.MixedCellTopology((10, 20, 30), 3)
    @test length(m_full) == 3
    @test collect(m_full) == [10, 20, 30]

    # invalid: len > MAX_K should throw on construction
    @test_throws ArgumentError ManifoldMeshes.MixedCellTopology((1, 2, 3), 4)
end

@testset "CellTypeStyle and PatchStyle for LatLonGrid" begin
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 90.0], lon_edges = [0.0, 360.0])
    @test CellTypeStyle(g) === IsUniform{4}()
    @test CellTypeStyle(LatLonGrid) === IsUniform{4}()
    @test PatchStyle(g) === NoPatch()
    @test PatchStyle(LatLonGrid) === NoPatch()
end
