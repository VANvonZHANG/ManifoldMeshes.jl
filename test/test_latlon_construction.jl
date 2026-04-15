@testset "LatLonGrid construction" begin
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 90.0, 180.0, 270.0, 360.0])

    @test g.nlat == 2
    @test g.nlon == 4
    @test g.R == 1.0

    # Global info
    @test manifold(g) isa Sphere
    @test num_cells(g) == 2 * 4
    @test num_nodes(g) == 3 * 5
    @test num_edges(g) == 3 * 4 + 2 * 4  # (nlat+1)*nlon + nlat*nlon

    # TopologyStyle
    @test TopologyStyle(g) === IsGrid()
    @test TopologyStyle(LatLonGrid) === IsGrid()
end

@testset "node_coordinates basics" begin
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 90.0, 180.0, 270.0, 360.0])
    # nlat=2, nlon=4 => nodes is (3, 5), nlon+1=5 columns per row

    # South pole: all nodes at lat=-90 map to (0, 0, -R)
    for ilon in 1:5
        n = node_coordinates(g, (1-1)*5 + ilon)
        @test n ≈ SVector(0.0, 0.0, -1.0) atol=1e-12
    end

    # North pole: all nodes at lat=90 map to (0, 0, R)
    for ilon in 1:5
        n = node_coordinates(g, (3-1)*5 + ilon)
        @test n ≈ SVector(0.0, 0.0, 1.0) atol=1e-12
    end

    # Equator at lon=0: node(2, 1) = lat_edges[2]=0, lon_edges[1]=0
    n = node_coordinates(g, (2-1)*5 + 1)
    @test n ≈ SVector(1.0, 0.0, 0.0) atol=1e-12

    # Equator at lon=90: node(2, 2) = lat_edges[2]=0, lon_edges[2]=90
    n = node_coordinates(g, (2-1)*5 + 2)
    @test n ≈ SVector(0.0, 1.0, 0.0) atol=1e-12

    # Equator at lon=180: node(2, 3) = lat_edges[2]=0, lon_edges[3]=180
    n = node_coordinates(g, (2-1)*5 + 3)
    @test n ≈ SVector(-1.0, 0.0, 0.0) atol=1e-12
end

@testset "constructor validation" begin
    # Missing -90 start
    @test_throws ArgumentError LatLonGrid(lat_edges=[-45.0, 0.0, 90.0], lon_edges=[0.0, 360.0])

    # Missing 90 end
    @test_throws ArgumentError LatLonGrid(lat_edges=[-90.0, 0.0, 45.0], lon_edges=[0.0, 360.0])

    # Not ascending
    @test_throws ArgumentError LatLonGrid(lat_edges=[-90.0, 90.0, 0.0], lon_edges=[0.0, 360.0])

    # Missing 0 start for lon
    @test_throws ArgumentError LatLonGrid(lat_edges=[-90.0, 90.0], lon_edges=[90.0, 360.0])

    # Missing 360 end for lon
    @test_throws ArgumentError LatLonGrid(lat_edges=[-90.0, 90.0], lon_edges=[0.0, 180.0])

    # Not ascending lon
    @test_throws ArgumentError LatLonGrid(lat_edges=[-90.0, 90.0], lon_edges=[0.0, 360.0, 180.0])

    # Need at least 2 edges each
    @test_throws ArgumentError LatLonGrid(lat_edges=[-90.0], lon_edges=[0.0, 360.0])
    @test_throws ArgumentError LatLonGrid(lat_edges=[-90.0, 90.0], lon_edges=[0.0])
end
