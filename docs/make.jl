using ManifoldMeshes
using Documenter

DocMeta.setdocmeta!(ManifoldMeshes, :DocTestSetup, :(using ManifoldMeshes); recursive = true)

makedocs(;
    modules = [ManifoldMeshes],
    sitename = "ManifoldMeshes.jl",
    format = Documenter.HTML(;
        canonical = "https://VANvonZHANG.github.io/ManifoldMeshes.jl",
        edit_link = "main",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Theory" => "theory.md",
        "Grid Types" => [
            "LatLonGrid" => "grids/latlon.md",
            "CubedSphereGrid" => "grids/cubed_sphere.md",
            "ReducedGaussianGrid" => "grids/reduced_gaussian.md",
            "HEALPixGrid" => "grids/healpix.md",
        ],
        "Visualization" => "visualization.md",
        "API Reference" => "api.md",
    ]
)

deploydocs(;
    repo = "github.com/VANvonZHANG/ManifoldMeshes.jl",
    devbranch = "main"
)
