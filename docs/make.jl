using PhyloSummaries
using Documenter

DocMeta.setdocmeta!(PhyloSummaries, :DocTestSetup, :(using PhyloSummaries); recursive=true)

makedocs(;
    modules=[PhyloSummaries],
    authors="Cecile Ane <cecileane@users.noreply.github.com>, Joshua Justison <jjustison@users.noreply.github.com>, and contributors",
    sitename="PhyloSummaries.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaPhylo.github.io/PhyloSummaries.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaPhylo/PhyloSummaries.jl",
    devbranch="main",
)
