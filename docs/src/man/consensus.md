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

## consensus tree of blob

Given a set of input phylogenetic networks (which may be trees or not),
we can obtain a consensus of their "tree of blobs".

For a given network, the tree of blob is obtained by contracting each
blob (2-edge connected component) in the network into a single node.
We get a tree, whose edges are exactly the "cut edges" of the network:
the edges that, when removed from the network, cut it into disconnected
components.

The example below uses a set of 10 networks from a file that comes with the
package, all of which have a single reticulation:

```@repl
netfile = joinpath(dirname(pathof(PhyloSummaries)), "..",
    "test","bootstrapnets_h1.nwk");
netsample = readmultinewick(netfile);
length(netsample) # 10 networks
netsample[1] # first network. taxa: t1 through t6
```

Let's look at the first and second, say:

```@example
R"svg"(figname("bootstrapnets_h1_12.svg"), width=7, height=3) # hide
R"layout"([1 2])       # hide
R"par"(mar=[0,0,1,0])  # hide
plot(netsample[1], showedgelength=true);
R"mtext"("net 1");     # hide
plot(netsample[2], showedgelength=true);
R"mtext"("net 2");
R"dev.off()"; # hide
nothing # hide
```

![nets 1 & 2, from bootstrapnets_h1](../assets/figures/bootstrapnets_h1_12.svg)

The function [`consensus_treeofblobs`](@ref) returns a `NamedTuple` with the following keys:

- `tob`: the consensus tree of blobs.
- `taxa`: the list of taxon names.
- `blob_table`: support for each blob partition in the consensus.
- `circorder_table` (ToB only): support for circular orderings within each blob.
- `hybrid_table`: support for hybrid clades within blobs.
- `bipartition_table`: support for cut-edge bipartitions not adjacent to
  an interesting blob (i.e. non-redundant bipartitions).

The tables are row-oriented (NamedTuples of vectors) and can be converted
to DataFrames.

```@repl
res = consensus_treeofblobs(netsample);
keys(res)
tree = res[:tob]
writenewick(tree, internallabel=false) # star tree! 1 big blob
```

We can inspect the support for the blobs and other features.
Here, we use the `DataFrames` package to display the tables.

```@repl
using DataFrames
res[:taxa]
blb_df = DataFrame(res[:blob_table], copycols=false) # use same columns in memory
bip_df = DataFrame(res[:bipartition_table], copycols=false)
hyb_df = DataFrame(res[:hybrid_table], copycols=false)
cir_df = DataFrame(res[:circorder_table], copycols=false)
```

Because the consensus is a star tree (a single blob, no cut edges),
`bip_df` is empty here.
The `blb_df` table shows the support for each blob partition.
The `hyb_df` table shows which taxa were hybrid in the input networks.
The `cir_df` table shows the support for circular orderings within blobs.

We can visualize these supports on the consensus tree.
The plot below shows:

- Support for each blob at the blob node (in blue).
- Support for its lowest hybrid descendants at each blob's edges (in black).

```@example
R"svg"(figname("consensus_tob_support.svg"), width=9, height=3) # hide
R"layout"([1 2 3]);
R"par"(mar=[0,0,0.5,0]);
plot(tree, shownodenumber=true, showedgenumber=true, tipoffset=0.05);
R"mtext"("node & edge numbers", side=3, line=-1);
plot(tree, nodelabelcolor="deepskyblue",
    nodelabeladj=1.1, edgelabeladj=[.5, -0.2],
    nodelabel=select(blb_df, [:node, :support_partition]),
    edgelabel=select(hyb_df, [:edge, :support_hybrid]));
R"mtext"("blob (blue) & hybrid (black) support", side=3, line=-1);
rotate!(netsample[8], -3);
rotate!(netsample[8], -6); # to de-tangle crossing edges in plot below
plot(netsample[8]);
R"mtext"("input net 8", side=3, line=-1);
R"dev.off()"; # hide
nothing # hide
```

![consensus tree of blobs with support values](../assets/figures/consensus_tob_support.svg)
