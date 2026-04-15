@testset "spherical triangle area" begin
    # Octant of unit sphere: vertices at (1,0,0), (0,1,0), (0,0,1)
    A = SVector(1.0, 0.0, 0.0)
    B = SVector(0.0, 1.0, 0.0)
    C = SVector(0.0, 0.0, 1.0)
    area = ManifoldMeshes._spherical_triangle_area(1.0, A, B, C)
    @test area ≈ π / 2 atol=1e-12  # octant = π/2 steradians
end

@testset "lobe test: total area = 4πR²" begin
    for R in [1.0, 6371.0]
        g = LatLonGrid(lat_edges=collect(-90.0:10.0:90.0),
                       lon_edges=collect(0.0:15.0:360.0), R=R)
        total = sum(i -> cell_volume(g, i), 1:num_cells(g))
        @test total ≈ 4 * π * R^2 rtol=1e-10
    end
end

@testset "symmetry: same latitude band = same volume" begin
    g = LatLonGrid(lat_edges=collect(-90.0:10.0:90.0),
                   lon_edges=collect(0.0:15.0:360.0))
    for ilat in 1:g.nlat
        id1 = ManifoldMeshes._cell_linear_index(g, ilat, 1)
        id2 = ManifoldMeshes._cell_linear_index(g, ilat, g.nlon)
        @test cell_volume(g, id1) ≈ cell_volume(g, id2) atol=1e-14
    end
end

@testset "geographic comparison near equator" begin
    R = 6371.0
    # Grid with 10° lat bands and 10° lon cells; band 10 is 0° to 10°N
    g = LatLonGrid(lat_edges=collect(-90.0:10.0:90.0),
                   lon_edges=collect(0.0:10.0:360.0), R=R)

    lat1, lat2 = deg2rad(0.0), deg2rad(10.0)
    Δlon = deg2rad(10.0)
    geographic_area = R^2 * (sin(lat2) - sin(lat1)) * Δlon

    # Cell in the 10th latitude band (0° to 10°), first longitude cell
    equatorial_cell = ManifoldMeshes._cell_linear_index(g, 10, 1)
    geodesic_area = cell_volume(g, equatorial_cell)
    @test geodesic_area ≈ geographic_area rtol=0.01
end

@testset "polar cell volume correctness" begin
    R = 1.0
    # Two bands: south hemisphere (-90 to 0), north hemisphere (0 to 90)
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 120.0, 240.0, 360.0], R=R)

    south_vol = cell_volume(g, ManifoldMeshes._cell_linear_index(g, 1, 1))
    @test south_vol > 0.0
    @test !isnan(south_vol)

    north_vol = cell_volume(g, ManifoldMeshes._cell_linear_index(g, 2, 1))
    @test north_vol ≈ south_vol atol=1e-14

    # South hemisphere total = 3 cells * south_vol (3 longitude sectors)
    south_total = sum(ilon -> cell_volume(g, ManifoldMeshes._cell_linear_index(g, 1, ilon)), 1:g.nlon)
    @test south_total ≈ 2 * π * R^2 atol=1e-12  # hemisphere = 2πR²

    # North hemisphere total
    north_total = sum(ilon -> cell_volume(g, ManifoldMeshes._cell_linear_index(g, 2, ilon)), 1:g.nlon)
    @test north_total ≈ 2 * π * R^2 atol=1e-12

    # Both hemispheres sum to full sphere
    @test south_total + north_total ≈ 4 * π * R^2 atol=1e-12
end

@testset "cell_centroid basics" begin
    R = 1.0
    g = LatLonGrid(lat_edges=[-90.0, 90.0], lon_edges=[0.0, 120.0, 240.0, 360.0], R=R)

    for ilon in 1:g.nlon
        c = cell_centroid(g, ManifoldMeshes._cell_linear_index(g, 1, ilon))
        @test abs(norm(c) - R) < 1e-10
        @test abs(c[3]) < 0.1  # should be near equator
    end
end

@testset "cell_centroid on sphere surface" begin
    g = LatLonGrid(lat_edges=collect(-90.0:10.0:90.0),
                   lon_edges=collect(0.0:15.0:360.0), R=1.0)
    for i in 1:num_cells(g)
        c = cell_centroid(g, i)
        @test abs(norm(c) - 1.0) < 1e-10
    end
end
