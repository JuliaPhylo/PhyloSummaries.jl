# installation

To install Julia, see [here](http://julialang.org/downloads/).

To read & write phylogenies (networks or trees, admixture graphs), you should
install [PhyloNetworks](https://github.com/juliaphylo/PhyloNetworks.jl),
which PhyloSummaries depends on. To install it, see
[here](https://juliaphylo.github.io/PhyloNetworks.jl/dev/man/installation/#Installation).

To visualize phylogenies, install [PhyloPlots](https://github.com/juliaphylo/PhyloPlots.jl).
We can do so in the Julia REPL for example:
enter package mode with `]`, and:

```
add PhyloPlots
```
Or in julian mode:

```julia
using Pkg
Pkg.add("PhyloPlots")
```
