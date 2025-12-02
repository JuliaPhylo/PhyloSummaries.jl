```@meta
ShareDefaultModule = true
```
```@setup
using PhyloSummaries
using PhyloPlots, RCall # not below because it would generate "precompiling..." output
figpath = joinpath("..", "assets", "figures"); mkpath(figpath)
figname(x) = joinpath(figpath, x)
```

# consensus phylogenies

## consensus of phylogenetic trees

Given a set of input phylogenies that are all trees, we can get their
greedy consensus or their majority-rule consensus with
[`consensustree`](@ref).

To give an example, we will use a set of trees from an example file
that comes with the package:

```@repl
using PhyloNetworks
inputfile = joinpath(dirname(pathof(PhyloSummaries)), "..","test","raxmltrees.tre");
treesample = readmultinewick(inputfile);
length(treesample) # 30 trees
treesample[1] # first tree in the list
```

To visualize trees and network, we use package
[PhyloPlots](https://github.com/juliaphylo/PhyloPlots.jl).

```julia
using PhyloPlots
using RCall            # to tweak our plot within R
```

```@example
R"svg"(figname("raxmltree_12.svg"), width=7, height=3) # hide
R"layout"([1 2])       # figure of 2 panels
R"par"(mar=[0,0,1,0])  # for smaller margins
plot(treesample[1], showedgelength=true);
R"mtext"("tree 1")     # add text annotation: title here
plot(treesample[2], showedgelength=true);
R"mtext"("tree 2")
R"dev.off()"; # hide
nothing # hide
```
![trees 1-2, from raxml.tre](../assets/figures/raxmltree_12.svg)

By default, we get the greedy consensus tree of our input trees,
considered as unrooted trees.

```@repl
con = consensustree(treesample)
writenewick(con, support=true)
```

To plot the consensus tree showing support values, we can first extract
the support values into a data frame, then use it to label edges.
Below, we multiple support values by 100 to get percentages.

```@repl
using DataFrames
esup = DataFrame(
    number = [e.number for e in con.edge if !isexternal(e)],
    support = [round(100 * e.y, digits=1) for e in con.edge if !isexternal(e)]
)
ecol = Dict(r[:number] => (r[:support] < 70 ? "orange2" : "black") for r in eachrow(esup))
```

```@example
R"svg"(figname("raxmltree_con.svg"), width=7, height=3) # hide
R"layout"([1 2])  # hide
R"par"(mar=[0,0,1,0])  # hide
plot(con, showedgenumber=true);
R"mtext"("edge numbers", line=0)
plot(con, edgelabel=esup, edgecolor=ecol);
R"mtext"("support: % input trees", line=0)
R"dev.off()"; # hide
nothing # hide
```
![majority rule consensus tree](../assets/figures/raxmltree_con.svg)

