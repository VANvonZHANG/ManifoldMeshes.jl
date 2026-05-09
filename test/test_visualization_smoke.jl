using ManifoldMeshes
using Manifolds
using Test

using CairoMakie

@testset "Visualization smoke tests" begin
    grid = LatLonGrid(lat_edges = collect(-90.0:30.0:90.0), lon_edges = collect(0.0:60.0:360.0))

    @testset "plot_mesh 3D wireframe" begin
        fig = plot_mesh(grid; show_nodes = false, show_edges = true,
            view = Symbol("3d"), figsize = (400, 400))
        @test fig isa Figure
        save("test_mesh_3d.png", fig)
        @test isfile("test_mesh_3d.png")
        rm("test_mesh_3d.png")
    end

    @testset "plot_mesh 2D equirectangular" begin
        fig = plot_mesh(grid; view = Symbol("2d"), figsize = (400, 400))
        @test fig isa Figure
        save("test_mesh_2d.png", fig)
        @test isfile("test_mesh_2d.png")
        rm("test_mesh_2d.png")
    end

    @testset "plot_mesh_filled 3D" begin
        fig = plot_mesh_filled(grid; show_edges = true, view = Symbol("3d"), figsize = (
            400, 400))
        @test fig isa Figure
        save("test_filled_3d.png", fig)
        @test isfile("test_filled_3d.png")
        rm("test_filled_3d.png")
    end

    @testset "plot_mesh_filled 2D" begin
        fig = plot_mesh_filled(grid; view = Symbol("2d"), figsize = (400, 400))
        @test fig isa Figure
        save("test_filled_2d.png", fig)
        @test isfile("test_filled_2d.png")
        rm("test_filled_2d.png")
    end

    @testset "invalid view raises error" begin
        @test_throws ErrorException plot_mesh(grid; view = :invalid)
    end

    @testset "color_by function" begin
        fig = plot_mesh_filled(grid; color_by = i -> (i % 2 == 0 ? :red : :blue),
            view = Symbol("3d"), figsize = (400, 400))
        @test fig isa Figure
    end
end
