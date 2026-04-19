@testset "edge_length at 45°N latitude" begin
    R = 1.0
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 45.0, 90.0],
        lon_edges = [0.0, 90.0, 180.0, 270.0, 360.0], R = R)

    # Meridional edge from 0°N to 45°N: angular length = π/4
    n_h = (g.nlat + 1) * g.nlon
    meridional = n_h + (2 - 1) * g.nlon + 1  # v(ilat=2, ilon=1): 0°→45°N
    @test edge_length(g, meridional) ≈ π / 4 atol=1e-8

    # Zonal edge at 45°N spanning 90°: geodesic (great circle) distance = π/3
    zonal_45 = (3 - 1) * g.nlon + 1  # h(ilat=3, ilon=1): at 45°N
    @test edge_length(g, zonal_45) ≈ π / 3 atol=1e-8
end

@testset "edge_midpoint degenerate polar edge" begin
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 90.0],
        lon_edges = [0.0, 90.0, 180.0, 270.0, 360.0])

    # South pole horizontal edge (edge_id=1): all nodes at south pole
    mp = edge_midpoint(g, 1)
    south_pole = SVector(0.0, 0.0, -1.0)
    @test mp ≈ south_pole atol=1e-12

    # North pole horizontal edge: all nodes at north pole
    n_h = (g.nlat + 1) * g.nlon
    mp_north = edge_midpoint(g, n_h)
    north_pole = SVector(0.0, 0.0, 1.0)
    @test mp_north ≈ north_pole atol=1e-12
end

@testset "fine grid (1° resolution) area conservation" begin
    R = 1.0
    g = LatLonGrid(lat_edges = collect(-90.0:1.0:90.0),
        lon_edges = collect(0.0:1.0:360.0), R = R)

    total = sum(i -> cell_volume(g, i), 1:num_cells(g))
    @test total ≈ 4 * π * R^2 rtol=1e-10
end

@testset "node_cells periodic boundary" begin
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 90.0],
        lon_edges = [0.0, 90.0, 180.0, 270.0, 360.0])

    # Node at equator, lon=0 (ilat=2, ilon=1)
    nid_0 = ManifoldMeshes._node_linear_index(g, 2, 1)
    cells_0 = sort(node_cells(g, nid_0))
    @test length(cells_0) == 4

    # Node at equator, lon=360 (ilat=2, ilon=5) — same physical point as lon=0
    nid_360 = ManifoldMeshes._node_linear_index(g, 2, 5)
    cells_360 = sort(node_cells(g, nid_360))
    @test length(cells_360) == 4

    # Both nodes should see the same cells (by linear index)
    @test cells_0 == cells_360
end

@testset "edge_outward_normal equatorial cells" begin
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 90.0],
        lon_edges = [0.0, 120.0, 240.0, 360.0])

    # Equatorial cell (ilat=2, ilon=1)
    eq_cell = ManifoldMeshes._cell_linear_index(g, 2, 1)
    cc = cell_centroid(g, eq_cell)

    for edge_id in cell_edges(g, eq_cell)
        result = edge_outward_normal(g, edge_id, eq_cell)
        bp, n = result.base_point, result.normal

        norm(n) < 1e-14 && continue

        # Orthogonal to radial direction (tangent space constraint)
        @test abs(dot(bp, n)) < 1e-10

        # Unit length
        @test abs(norm(n) - 1.0) < 1e-10

        # Points outward from cell center
        diff = cc - bp
        proj_diff = diff - dot(diff, bp) * bp
        @test dot(n, proj_diff) < 0
    end
end
