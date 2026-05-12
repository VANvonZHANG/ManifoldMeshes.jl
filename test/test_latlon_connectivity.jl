@testset "cell_nodes consistency" begin
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 90.0], lon_edges = [
        0.0, 90.0, 180.0, 270.0, 360.0])

    # Spot check one cell (SW corner)
    sn = cell_nodes(g, 1)
    @test length(sn) == 4
    @test node_coordinates(g, sn[1]) ≈ SVector(0.0, 0.0, -1.0) atol=1e-12
    # Return type
    @test cell_nodes(g, 1) isa NTuple{4, Int}
end

@testset "cell_cells periodicity" begin
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 90.0], lon_edges = [
        0.0, 90.0, 180.0, 270.0, 360.0])

    # Cell at (ilat=1, ilon=1): west neighbor is (1, nlon=4)
    _, _, west, _ = cell_cells(g, 1)
    @test west == ManifoldMeshes._cell_linear_index(g, 1, 4)

    # Cell at (ilat=1, ilon=4): east neighbor is (1, 1)
    _, _, _, east = cell_cells(g, ManifoldMeshes._cell_linear_index(g, 1, 4))
    @test east == ManifoldMeshes._cell_linear_index(g, 1, 1)
end

@testset "cell_cells sentinel values" begin
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 90.0], lon_edges = [0.0, 360.0])

    south, north, west, east = cell_cells(g, 1)
    @test south == 0
    @test north == 2
    @test west == 1  # periodic (only 1 column)
    @test east == 1  # periodic

    south, north, west, east = cell_cells(g, 2)
    @test south == 1
    @test north == 0

    @test cell_cells(g, 1) isa NTuple{4, Int}
end

@testset "node_cells interior vs pole" begin
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 90.0], lon_edges = [
        0.0, 90.0, 180.0, 270.0, 360.0])

    # Interior node (ilat=2, ilon=1): equator, lon=0 — 4 adjacent cells
    nc = node_cells(g, ManifoldMeshes._node_linear_index(g, 2, 1))
    @test length(nc) == 4

    # South pole node (ilat=1, ilon=1): 2 adjacent cells
    nc = node_cells(g, ManifoldMeshes._node_linear_index(g, 1, 1))
    @test length(nc) == 2

    # North pole node (ilat=3, ilon=1): 2 adjacent cells
    nc = node_cells(g, ManifoldMeshes._node_linear_index(g, 3, 1))
    @test length(nc) == 2
end

@testset "cell_edges consistency with cell_nodes" begin
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 90.0], lon_edges = [
        0.0, 90.0, 180.0, 270.0, 360.0])

    # Spot check one cell
    cell_id = 1
    edges = cell_edges(g, cell_id)
    @test edges isa NTuple{4, Int}

    sn = cell_nodes(g, cell_id)
    n1, n2 = ManifoldMeshes._edge_endpoints(g, edges[1])  # south edge
    @test n1 == node_coordinates(g, sn[1])  # SW
    @test n2 == node_coordinates(g, sn[2])  # SE
end

@testset "cell_edges ordering on interior cell" begin
    g = LatLonGrid(lat_edges = collect(-90.0:30.0:90.0),
        lon_edges = collect(0.0:60.0:360.0))
    # Interior cell (ilat=2, ilon=2) — not at pole or periodic boundary
    cell_id = ManifoldMeshes._cell_linear_index(g, 2, 2)
    edges = cell_edges(g, cell_id)
    sn = cell_nodes(g, cell_id)

    # South edge: SW → SE
    n1, n2 = ManifoldMeshes._edge_endpoints(g, edges[1])
    @test n1 == node_coordinates(g, sn[1])
    @test n2 == node_coordinates(g, sn[2])

    # North edge: NW → NE
    n1, n2 = ManifoldMeshes._edge_endpoints(g, edges[2])
    @test n1 == node_coordinates(g, sn[4])
    @test n2 == node_coordinates(g, sn[3])

    # West edge: SW → NW
    n1, n2 = ManifoldMeshes._edge_endpoints(g, edges[3])
    @test n1 == node_coordinates(g, sn[1])
    @test n2 == node_coordinates(g, sn[4])

    # East edge: SE → NE
    n1, n2 = ManifoldMeshes._edge_endpoints(g, edges[4])
    @test n1 == node_coordinates(g, sn[2])
    @test n2 == node_coordinates(g, sn[3])
end

@testset "cell_edges periodicity" begin
    g = LatLonGrid(lat_edges = [-90.0, 0.0, 90.0], lon_edges = [
        0.0, 90.0, 180.0, 270.0, 360.0])

    # Cell (1,1) and cell (1,4) share the western edge of cell (1,1) / eastern edge of cell (1,4)
    _, _, w1, _ = cell_edges(g, ManifoldMeshes._cell_linear_index(g, 1, 1))
    _, _, _, e4 = cell_edges(g, ManifoldMeshes._cell_linear_index(g, 1, 4))
    @test w1 == e4  # same edge, periodic wrap
end
