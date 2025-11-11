using PhyloSummaries
using Documenter

using DocumenterInterLinks
links = InterLinks(
    "PhyloNetworks" => "https://juliaphylo.github.io/PhyloNetworks.jl/stable/objects.inv",
)

# loading PhyloNetworks by default in all docstring examples
DocMeta.setdocmeta!(
    PhyloSummaries,
    :DocTestSetup,
    :(using PhyloNetworks; using PhyloSummaries);
    recursive=true)

makedocs(;
    modules=[PhyloSummaries],
    authors="Cecile Ane <cecileane@users.noreply.github.com>, Ayush Sharma <ayusharma17@users.noreply.github.com>, Joshua Justison <jjustison@users.noreply.github.com>, and contributors",
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
    push_preview = true,
)
