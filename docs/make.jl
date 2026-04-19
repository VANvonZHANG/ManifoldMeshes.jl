using ManifoldMeshes
using Documenter

DocMeta.setdocmeta!(ManifoldMeshes, :DocTestSetup, :(using ManifoldMeshes); recursive=true)

makedocs(;
    modules=[ManifoldMeshes],
    sitename="ManifoldMeshes.jl",
    format=Documenter.HTML(;
        canonical="https://VANvonZHANG.github.io/ManifoldMeshes.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/VANvonZHANG/ManifoldMeshes.jl",
    devbranch="main",
)
