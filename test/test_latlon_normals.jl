using LinearAlgebra

@testset "edge_length equatorial" begin
    R = 1.0
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 90.0, 180.0, 270.0, 360.0], R=R)

    # Vertical edge at south pole to equator: length = pi/2
    n_h = (g.nlat + 1) * g.nlon
    south_to_eq = n_h + 1  # v(1,1)
    @test edge_length(g, south_to_eq) ≈ π/2 atol=1e-12

    # Horizontal edge at equator spanning 90 degrees: length = pi/2
    eq_edge = (2-1)*g.nlon + 1  # h(2,1)
    @test edge_length(g, eq_edge) ≈ π/2 atol=1e-12
end

@testset "edge_midpoint equidistant" begin
    g = LatLonGrid(lat_edges=[-90.0, 90.0], lon_edges=[0.0, 180.0, 360.0], R=1.0)

    for edge_id in 1:num_edges(g)
        mp = edge_midpoint(g, edge_id)
        n1, n2 = ManifoldMeshes._edge_endpoints(g, edge_id)
        d1 = Manifolds.distance(g.manifold, mp, n1)
        d2 = Manifolds.distance(g.manifold, mp, n2)
        @test d1 ≈ d2 atol=1e-12
        @test abs(norm(mp) - 1.0) < 1e-10
    end
end

@testset "edge_outward_normal orthogonality" begin
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 90.0, 180.0, 270.0, 360.0])

    for cell_id in 1:num_cells(g)
        edges = cell_edges(g, cell_id)
        cc = cell_centroid(g, cell_id)
        for edge_id in edges
            result = edge_outward_normal(g, edge_id, cell_id)
            bp, n = result.base_point, result.normal

            # Skip degenerate (polar collapse) edges
            norm(n) < 1e-14 && continue

            # Normal should be in tangent space: dot(bp, n) ~= 0
            @test abs(dot(bp, n)) < 1e-10

            # Normal should be unit length
            @test abs(norm(n) - 1.0) < 1e-10

            # Normal should point outward: dot(normal, centroid-base_point tangent projection) < 0
            diff = cc - bp
            proj_diff = diff - dot(diff, bp) * bp  # project to tangent space at bp
            @test dot(n, proj_diff) < 0
        end
    end
end

@testset "polar collapse: degenerate edge normal" begin
    g = LatLonGrid(lat_edges=[-90.0, 0.0, 90.0], lon_edges=[0.0, 90.0, 180.0, 270.0, 360.0])

    # South pole horizontal edge: all nodes at south pole, zero-length edge
    pole_edge = 1  # h(1,1)
    result = edge_outward_normal(g, pole_edge, 1)
    @test result.normal == zero(SVector{3,Float64})
    @test !any(isnan, result.normal)

    # North pole horizontal edge
    n_h = (g.nlat + 1) * g.nlon
    north_pole_edge = n_h
    result = edge_outward_normal(g, north_pole_edge, ManifoldMeshes._cell_linear_index(g, 2, 1))
    @test result.normal == zero(SVector{3,Float64})
    @test !any(isnan, result.normal)
end

@testset "boundary markers (full sphere = no boundary)" begin
    g = LatLonGrid(lat_edges=[-90.0, 90.0], lon_edges=[0.0, 360.0])
    @test boundary_nodes(g, :any) == Int[]
    @test boundary_edges(g, :any) == Int[]
    @test boundary_nodes(g, 0) == Int[]
    @test boundary_edges(g, 0) == Int[]
end
