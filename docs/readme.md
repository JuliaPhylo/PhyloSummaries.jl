## to build the documentation locally

Run the commands below to:
- test the `jldoctest` blocks of examples in the docstrings
- create or update a `build/` directory with html files

```shell
cd docs
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.develop(PackageSpec(path=pwd()))'
julia --project --color=yes make.jl
```

or interactively in `docs/`:

```shell
pkg> activate .
pkg> instantiate
pkg> dev ~/.julia/dev/PhyloSummaries
julia> include("make.jl")
```

Then open `src/build/index.html`.
