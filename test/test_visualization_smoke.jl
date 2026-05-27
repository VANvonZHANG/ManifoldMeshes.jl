using ManifoldMeshes
using Manifolds
using Test

using CairoMakie

@testset "Visualization smoke tests" begin
    grids = [
        ("LatLonGrid", LatLonGrid(lat_edges = collect(-90.0:30.0:90.0), lon_edges = collect(0.0:60.0:360.0))),
        ("CubedSphereGrid", CubedSphereGrid(n = 4)),
        ("ReducedGaussianGrid", ReducedGaussianGrid(nlat = 8)),
        ("HEALPixGrid", HEALPixGrid(nside = 2)),
    ]

    for (name, grid) in grids
        @testset "plot_mesh 3D wireframe on $name" begin
            fig = plot_mesh(grid; show_nodes = false, show_edges = true, figsize = (400, 400))
            @test fig isa Figure
            save("test_mesh_$(name).png", fig)
            @test isfile("test_mesh_$(name).png")
            rm("test_mesh_$(name).png")
        end

        @testset "plot_mesh_filled 3D on $name" begin
            fig = plot_mesh_filled(grid; show_edges = true, figsize = (400, 400))
            @test fig isa Figure
            save("test_filled_$(name).png", fig)
            @test isfile("test_filled_$(name).png")
            rm("test_filled_$(name).png")
        end
    end

    @testset "color_by function" begin
        grid = LatLonGrid(lat_edges = collect(-90.0:30.0:90.0), lon_edges = collect(0.0:60.0:360.0))
        fig = plot_mesh_filled(grid; color_by = i -> (i % 2 == 0 ? :red : :blue), figsize = (400, 400))
        @test fig isa Figure
    end
end
