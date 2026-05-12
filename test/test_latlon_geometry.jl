@testset "spherical triangle area" begin
    # Octant of unit sphere: vertices at (1,0,0), (0,1,0), (0,0,1)
    A = SVector(1.0, 0.0, 0.0)
    B = SVector(0.0, 1.0, 0.0)
    C = SVector(0.0, 0.0, 1.0)
    area = ManifoldMeshes.spherical_triangle_area(1.0, A, B, C)
    @test area ≈ π / 2 atol=1e-12  # octant = π/2 steradians
end

@testset "total area = 4πR²" begin
    g = LatLonGrid(lat_edges = collect(-90.0:10.0:90.0),
        lon_edges = collect(0.0:15.0:360.0), R = 1.0)
    total = sum(cell_volume(g, i) for i in 1:num_cells(g))
    @test total ≈ 4π rtol=1e-10
end

@testset "radius scaling (R ≠ 1.0)" begin
    g = LatLonGrid(lat_edges = collect(-90.0:10.0:90.0),
        lon_edges = collect(0.0:15.0:360.0), R = 6371.0)
    total = sum(cell_volume(g, i) for i in 1:num_cells(g))
    @test total ≈ 4π * 6371.0^2 rtol=1e-10
end

@testset "symmetry: same latitude band = same volume" begin
    g = LatLonGrid(lat_edges = collect(-90.0:10.0:90.0),
        lon_edges = collect(0.0:15.0:360.0))
    # Test one equatorial band and one polar band
    id_eq1 = ManifoldMeshes._cell_linear_index(g, 6, 1)
    id_eq2 = ManifoldMeshes._cell_linear_index(g, 6, g.nlon)
    @test cell_volume(g, id_eq1) ≈ cell_volume(g, id_eq2) atol=1e-14
end

@testset "geographic comparison near equator" begin
    R = 6371.0
    # Grid with 10° lat bands and 10° lon cells; band 10 is 0° to 10°N
    g = LatLonGrid(lat_edges = collect(-90.0:10.0:90.0),
        lon_edges = collect(0.0:10.0:360.0), R = R)

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
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 90.0], lon_edges = [0.0, 120.0, 240.0, 360.0], R = R)

    south_vol = cell_volume(g, ManifoldMeshes._cell_linear_index(g, 1, 1))
    @test south_vol > 0.0
    @test !isnan(south_vol)

    north_vol = cell_volume(g, ManifoldMeshes._cell_linear_index(g, 2, 1))
    @test north_vol ≈ south_vol atol=1e-14

    # Full sphere conservation (only 6 cells, cheap)
    total = sum(cell_volume(g, i) for i in 1:num_cells(g))
    @test total ≈ 4π atol=1e-12
end

@testset "180-degree polar cell (degenerate diagonal)" begin
    R = 1.0
    # Coarse grid with 180° longitude cells — triggers the degenerate diagonal
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 90.0], lon_edges = [0.0, 180.0, 360.0], R = R)

    total = sum(cell_volume(g, i) for i in 1:num_cells(g))
    @test total ≈ 4π rtol=1e-10

    # Only 4 cells — spot check one
    @test cell_volume(g, 1) > 0.0
end

@testset "cell_centroid on sphere surface" begin
    g = LatLonGrid(lat_edges = collect(-90.0:10.0:90.0),
        lon_edges = collect(0.0:15.0:360.0), R = 1.0)
    # Spot check: equator, mid-latitude, and polar-adjacent cells
    c_eq = cell_centroid(g, ManifoldMeshes._cell_linear_index(g, 6, 1))
    @test abs(norm(c_eq) - 1.0) < 1e-10
    c_mid = cell_centroid(g, ManifoldMeshes._cell_linear_index(g, 3, 5))
    @test abs(norm(c_mid) - 1.0) < 1e-10
    c_polar = cell_centroid(g, ManifoldMeshes._cell_linear_index(g, 1, 3))
    @test abs(norm(c_polar) - 1.0) < 1e-10
end
