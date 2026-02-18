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
ecol = Dict(r[:number] => (r[:support] < 70 ? "orange2" : "black")
            for r in eachrow(esup))
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

- `tob`: the consensus tree of blobs
- `taxa`: list of taxon names, alphabetically, to connect with
  taxon numbers in some columns of the tables below
- `blob_table`: support for each blob partition in the consensus.
- `circorder_table`: support for circular orderings of the
  taxon blocks within each blob
- `hybrid_table`: support for hybrid clades within blobs
- `bipartition_table`: support for cut-edge bipartitions not adjacent to
  an interesting blob (i.e. non-redundant bipartitions).

The tables are row-oriented (NamedTuples of vectors) and can be converted
to data frames.

```@repl
res = consensus_treeofblobs(netsample);
keys(res)
tree = res[:tob]
writenewick(tree, internallabel=false) # star tree! 1 big blob
writenewick(tree) # internal node names: indicate their numbers used later
```

Here our consensus tree of blob is a star!
The root is a polytomy with 6 edges: one to each taxon.
This root node represents a blob that is present in many input networks.
It was assigned the name `_7_blob1`, so later this node will be referred
as "node 7" and stands for "blob 1".

The other parts of the results are tables that
let us inspect the support for the blob(s) and other features.
Here, we use the [`DataFrames`](https://dataframes.juliadata.org/stable/)
package to display the tables.

```@repl
using DataFrames
res[:taxa]
blb_df = DataFrame(res[:blob_table], copycols=false); # re-use in memory
show(blb_df, allcols=true)
```

This `blb_df` table describes each blob partition, and shows its support.
Here our consensus tree of blobs has a single blob, so this table
has a single row.
- When the edges and nodes of this blob are removed, the full taxon set
  is disconnected into multiple taxon blocks: a *partition* of the taxa.
  This partition is described by a list of taxa in each block,
  with blocks separated by `|`. So here we see that our partition has
  6 parts, or *taxon blocks*, and each taxon is alone in its own block.
- The number of taxon blocks is the **degree** of the blob, given in
  the "degree" column.
- The "node" column given the number of the node that represents the blob, in
  the tree of blobs. It's a node with ≥ 4 edges (a polytomy), "degree" edges
  actually: here 6 edges. Here, it's the root. This node number is indicated
  in the name assigned to internal node (here 7: in `_7_blob1`, see above).
- The `support_partition` column gives the proportion of input networks
  that have a blob with this particular partition: when the blob is removed
  from the input network, the taxa are divided into these same blocks.
  Here, 60% of our networks have this partition --and 40% don't.
  Scrolling back up, we see that "net 1" does have a blob with this partition,
  and "net 2" does not. This second network "net 2" has `t1,t2` together in
  the same block, when we delete its blob.
- The partition is given twice: with taxon numbers in column "partition_num",
  and with taxon names in column "partition".
  With longer taxon names than in this simple example, the full partition
  is often long and not shown in full. To see more with all taxon names,
  we can increase the `truncate` number. For example:

```@repl
show(select(blb_df, [:support_partition,:partition]), truncate=64)
```

Next we can look the table describing bipartitions: edges in the
consensus tree of blobs that, when deleted, separate the taxon into
*two* parts. Since many edges are already in the tree because they
are connected to a blob, we only look at the support of edges as being
non-redundant with any blob, that is, *additional* support: the proportion
of input networks that have this edge *not* connected to any blob.

Here we don't have any such edges because our consensus is a
star tree (a single blob, no cut edges).
So the corresponding table `bip_df` is empty:

```@repl
bip_df = DataFrame(res[:bipartition_table], copycols=false)
```

More interestingly then, we can look at which taxon blocks are
coming off of a (any) blob from a hybrid node:

```@repl
hyb_df = DataFrame(res[:hybrid_table], copycols=false)
```

This `hyb_df` table shows which taxon blocks are of hybrid origin
in the input networks, and in which proportion of networks.
A taxon block is of "hybrid origin" if it is below a hybrid node
that "exits" a blob: a hybrid with a child edge outside the blob
that the hybrid is in, and whose descendant clade is exactly this
taxon block.

In our example:
- The most frequent hybrid is t6: in 50% of input networks.
  For example, t6 is of hybrid origin in "net 1" above.
  Note that it may be of hybrid origin in networks that do not
  have the blob in our consensus tree of blobs.
- The other rows list the other taxon blocks of hybrid origin
  in some input networks: t3, t4 (each in 20% of networks) and t1 (10%).
  For example, t3 is of hybrid origin in "net 2", although below
  a blob that partitions the taxa differently than the blob
  retained in the consensus.

The table about circular orders is about blobs of level 1 in input networks,
for which the taxon blocks are arranged around the blob like around a circle.
This "circular order" can be read clock-wise or counter-clockwise, and
starting from any taxon block.
The table below lists all the circular orders found in the input network,
for the blob partitions retained in the consensus.

```@repl
cir_df = DataFrame(res[:circorder_table], copycols=false);
show(cir_df, allcols=true)
```

Here, this table has a single row: the same circular order was found in
all input networks that have the blob in our consensus
(60% of input networks from the blob information in `blb_df` above). This
circular order is indicated in the columns `partition_num` (taxon numbers)
and `partition` (taxon names). Turning around the cycle, and starting
arbitrarily at t2, our networks have this order:
t2 → t4 → t3 → t5 → t6 → t1 (then back to t2).

Below, we visualize some of these support values on the consensus
tree of blobs. We show:

- in blue, at nodes: the support for the blob partition at each blob node
- in black, along edges adjacent to a blob: the support for each taxon block
  to be of hybrid origin (in input networks that have this blob or not).

To annotate the appropriate nodes and edges of the consensus,
we use the node and edge numbers in the table. To clarify, we first plot
the consensus with nodes and edges annotated with their numbers.
And for fun, we plot the 8th intput network: which has the blob in
the consensus, but with t1 of hybrid origin, not t6 like in "net 1".

```@example
R"svg"(figname("consensus_tob_support.svg"), width=9, height=3) # hide
R"layout"([1 2 3]);
R"par"(mar=[0,0,0.5,0]);
plot(tree, shownodenumber=true, showedgenumber=true, tipoffset=0.05);
R"mtext"("node & edge numbers", side=3, line=-1);
plot(tree, nodelabelcolor="orangered", edgelabelcolor="deepskyblue",
    nodelabeladj=1.1, edgelabeladj=[.5, -0.2],
    nodelabel=select(blb_df, [:node, :support_partition]),
    edgelabel=select(hyb_df, [:edge, :support_hybrid]));
R"mtext"("blob (red) & hybrid (blue) support", side=3, line=-1);
rotate!(netsample[8], -3);
rotate!(netsample[8], -6); # to de-tangle crossing edges in plot below
plot(netsample[8]);
R"mtext"("input net 8", side=3, line=-1);
R"dev.off()"; # hide
nothing # hide
```

![consensus tree of blobs with support values](../assets/figures/consensus_tob_support.svg)

## consensus level-1 network

The function [`consensus_level1network`](@ref) first calculates
the consensus tree of blobs, then expands each blob node in that tree,
choosing the best-supported circular ordering and hybrid clade for each blob.
This function only takes level-1 networks as input: in which each blob
is a cycle with a circular order and a single hybrid node.
The result is a level-1 phylogenetic network.

Using the same set of input networks:

```@repl
res_l1 = consensus_level1network(netsample);
keys(res_l1)
con = res_l1[:net] # node labels match node numbers used later
writenewick(con, internallabel=false)
```

The plot below shows our consensus level-1 network, with node and edges
annotated by their numbers. We used the first plot to decide at which
nodes to rotate edges, to untangle the crossing lines.

```@example
R"svg"(figname("consensus_l1net_names.svg"), width=7, height=3) # hide
R"layout"([1 2])       # hide
R"par"(mar=[0,0,0.5,0]);
plot(con, shownodenumber=true, showedgenumber=true, tipoffset=0.05);
rotate!(con, 12); # rotating edges below node 12 will un-cross lines
R"mtext"("node & edge numbers", side=3, line=-1);
plot(con, shownodelabel=true);
R"mtext"("node names: match node numbers", side=3, line=-1);
R"dev.off()"; # hide
nothing # hide
```

![consensus level-1 network with edge/node numbers and names](../assets/figures/consensus_l1net_names.svg)

Like for the consensus tree of blobs, the rest of the result gives
us support values for various features of our consensus.

The `blob_table` for a level-1 network consensus includes additional columns
compared to the tree-of-blobs version: `hybrid` (the hybrid node number),
`support_circorder` (support for the chosen circular ordering) and
`hybrid_cluster` (the taxon block chosen as the hybrid clade):

```@repl
blb_l1 = DataFrame(res_l1[:blob_table], copycols=false)
names(blb_l1) # column names: including those not shown above
select(blb_l1, r"hybrid") # columns with 'hybrid' in col name
show(select(blb_l1, r"partition"), truncate=64)
```

We know from the tree-of-blobs consensus that other taxon blocks (like t1)
are of hybrid origin in some input networks.
The `blb_l1` table only lists t6 as hybrid, because it's the hybrid used
to build the level-1 consensus network, and it has 1 row per blob.
We see that t6 is the block below the blob's hybrid node in the consensus,
and it is of hybrid origin (perhaps in different types of blobs)
in 60% of input networks.

But we get the information that t1 is an alternative from the table
about hybrid clades, which can have multiple rows per blob,
with the same columns as from the tree-of-blobs consensus above:

```@repl
hyb_l1 = DataFrame(res_l1[:hybrid_table], copycols=false)
```

It's easier to see support values on a plot. Below, no edge has any
additional support for being non-redundant with a blob (because our
blob is large: no internal edge is outside the blob).

```@example
R"svg"(figname("consensus_l1net.svg"), width=7, height=3) # hide
R"layout"([1 2])       # hide
R"par"(mar=[0,0,0.5,0]); # hide
plot(con, shownodenumber=true, showedgenumber=true, tipoffset=0.05);
R"mtext"("node & edge numbers", side=3, line=-.5);
plot(con, nodelabelcolor="orangered", edgelabelcolor="deepskyblue",
    nodelabeladj=-0.1, edgelabeladj=[.5, -0.2],
    nodelabel=select(blb_l1, [:node, :support_partition]),
    edgelabel=select(hyb_l1, [:edge, :support_hybrid]));
R"mtext"("blob (red) & hybrid (blue) support", side=3, line=-.5);
R"dev.off()"; # hide
nothing # hide
```

![consensus level-1 network with support values](../assets/figures/consensus_l1net.svg)

### saving the results

To save the consensus network and its data tables to files, use
[`consensus_level1network_save`](@ref). The following will create
4 files: `level1consensus_net.nwk` with the consensus network in newick format,
and the various tables in csv format:
`level1consensus_blob.csv`, `level1consensus_hybrid.csv` and
`level1consensus_bipartition.csv`.

```@repl
consensus_level1network_save(res_l1, "level1consensus")
```

We can re-read these files later, like this:

```@example
con2 = readnewick("level1consensus_net.nwk");
writenewick(con2) # internal node names give node numbers used in tables
```
```@example
using CSV
blb2 = CSV.read("level1consensus_blob.csv", DataFrame)
hyb2  = CSV.read("level1consensus_hybrid.csv", DataFrame); # no edge column
bip2  = CSV.read("level1consensus_bipartition.csv", DataFrame); # no edge
nothing # hide
```

However, care should be taken after re-reading the network from the newick
file to recover the original node and edge numbers, because `readnewick` may set
them differently. We want to recover internal node numbers in the network
that match with those saved in the csv files.
To get the original node numbers, which are part of the node labels,
use [`resetnodenumbers_fromnames!`](@ref) like below.

```@repl
printnodes(con2) # node numbers are in the first column
resetnodenumbers_fromnames!(con2)
printnodes(con2) # new node numbers: now match names. everything else the same
```

Now, we can re-do our plot, but the original edge numbers are not recovered.
To map hybrid support, or support for edges as coming off a hybrid node,
we need to get the current number of each edge based on its nodes.

```@repl
hyb2 # no 'edge' column, because its old numbers are different after readnewick
hyb2.edge = edgenumbers_fromnodenumbers(hyb2, con2);
hyb2 # new edge column: with number that match those in "con2" network
```

Edge support, as non-redundant bipartitions, can be extracted from the network
itself: these support values were saved in the newick file.
Or it can be recovered from the bipartition table like we did for the
hybrid table above:

```@example
bip2.edge = edgenumbers_fromnodenumbers(bip2, con2)
nothing # hide
```

Here this is not interesting because no edge had support as
a non-redundant bipartition, so this table is empty.

```@example
R"svg"(figname("consensus_l1net2.svg"), width=7, height=3) # hide
R"layout"([1 2])       # hide
R"par"(mar=[0,0,0.5,0]); # hide
plot(con2, shownodenumber=true, showedgenumber=true, tipoffset=0.05);
R"mtext"("node & edge numbers", side=3, line=-.5);
R"mtext"("nodes: as before. edges: different", side=1, line=-1);
plot(con2, nodelabelcolor="orangered", edgelabelcolor="deepskyblue",
    nodelabeladj=-0.1, edgelabeladj=[.5, -0.2],
    nodelabel = select(blb2, [:node, :support_partition]),
    edgelabel = select(hyb2, [:edge, :support_hybrid])
);
R"mtext"("blob (red) & hybrid (blue) support", side=3, line=-.5);
R"dev.off()"; # hide
nothing # hide
```

![consensus level-1 network with support values, after re-reading results from file](../assets/figures/consensus_l1net2.svg)


```@eval
for f in "level1consensus_" .* ["net.nwk","blob.csv","hybrid.csv","bipartition.csv"]
    rm(f)
end # clean-up
nothing
```

### outgroup

We can specify one outgroup to root the consensus network. This forces the root
to be along the edge between the outgroup taxon and the rest of the network,
and prevents any reticulation from the ingroup into this outgroup.

Here, for example, if we indicate that t6 is the outgroup, then the second most
frequent hybrid block is chosen instead of t6. We saw that t3 and t4 are tied
next (in section [consensus tree of blob](@ref)).
One of the 2 will be chosen arbitrarily, here t4:

```@example
res_t6out = consensus_level1network(netsample; outgroup="t6", suppressinfo=true);
writenewick(res_t6out[:net], internallabel=false) # t4 below the hybrid H10
```

```@example
R"svg"(figname("consensus_l1net_t6out.svg"), width=7, height=3) # hide
R"layout"([1 2])       # hide
R"par"(mar=[0,0,0.5,0]); # hide
plot(res_t6out[:net], shownodenumber=true, tipoffset=0.05);
R"mtext"("node numbers", side=3, line=-.5);
R"mtext"("to see around which ones to rotate edges", side=1, line=-1); #TODO
for i in [13,9] rotate!(res_t6out[:net], i); end
plot(res_t6out[:net], tipoffset=0.05, nodelabeladj=-0.1, edgelabeladj=[.5,-0.2],
    nodelabelcolor="orangered", edgelabelcolor="deepskyblue",
    nodelabel = DataFrame(res_t6out[:blob_table][[:node, :support_partition]]),
    edgelabel = DataFrame(res_t6out[:hybrid_table][[:edge, :support_hybrid]])
);
R"mtext"("blob (red) & hybrid (blue) support", side=3, line=-.5);
R"mtext"("we see t6 hybrid in 50% input nets", side=1, line=-1);
R"dev.off()"; # hide
nothing # hide
```

![consensus level-1 network, required to be rooted with t6 as outgroup](../assets/figures/consensus_l1net_t6out.svg)

### non-redundant bipartitions

We use a different example for this: again a small set of small input networks,
so it's easy to look at them. They are all of level 1.

```@example
netfile = joinpath(dirname(pathof(PhyloSummaries)), "..",
    "test","level1_7taxa_abc.nwk");
netsample = readmultinewick(netfile);
length(netsample) # 5 networks: very small example
```

We first get the consensus tree of blobs, to check if there are multiple
circular orders, and if the top ones may be tied (we don't get this info
from the level-1 consensus):

```@repl
tob = consensus_treeofblobs(netsample, suppressinfo=true);
DataFrame(tob[:circorder_table]) # no ties: only 1 order
```

We see only 1 circular order, meaning that the input networks with the
top-ranking blob all share the same circular order.
In particular, there's no tie (good to know).

```@repl
con = consensus_level1network(netsample, suppressinfo=true);
writenewick(con[:net], support=true) # edges support: as non-redundant biparts
bip_df = DataFrame(con[:bipartition_table])
```

Let's plot our level-1 consensus: with support values on edges that show
support for hybrid clades (middle) or for
bipartitions to be non-redundant with blobs (right):

```@example
net = con[:net]
blb_df = DataFrame(con[:blob_table])
hyb_df = DataFrame(con[:hybrid_table])
bip_df = DataFrame(con[:bipartition_table])
R"svg"(figname("consensus_7taxa_abc.svg"), width=9, height=3) # hide
R"layout"([1 2 3]); # hide
R"par"(mar=[0,0,0.5,0]); # hide
plot(net, shownodenumber=true, tipoffset=0.05);
R"mtext"("node numbers: to decide rotations", side=3, line=-1);
for i in [14,13] rotate!(net, i); end
plot(net, nodelabeladj=1.1, edgelabeladj=[.5,-0.2],
    nodelabelcolor="orangered", edgelabelcolor="deepskyblue",
    nodelabel=select(blb_df, [:node, :support_partition]),
    edgelabel=select(hyb_df, [:edge, :support_hybrid]));
R"mtext"("blob (red) & hybrid (blue) support", side=3, line=-1);
plot(net, nodelabeladj=1.1, edgelabeladj=[.5,-0.2], nodelabelcolor="orangered",
    nodelabel=select(blb_df, [:node, :support_partition]),
    edgelabel=select(bip_df, [:edge, :support_nonredundant]));
R"mtext"("support for edges as non-redundant (black)", side=3, line=-1);
R"dev.off()"; # hide
nothing # hide
```

![consensus level-1 network with support values, 5 input networks on 7 taxa](../assets/figures/consensus_7taxa_abc.svg)

- (b1,c1,c2) has 20% support as a non-redundant bipartition.
  We also know that this clade is in the 60% input networks that have the blob
  (and redundant with the blob in these networks),
  so it's in at least 80% of the input networks.
- (c1,c2) has 60% support as a non-redundant bipartition. This clade may
  be present in more than 60% networks, if it's connected to a blob in these
  other networks, but we don't get to see this here.

### support for hybrid sisters

There may be high support for a clade to be of hybrid origin, in which case
this clade has 2 sisters (in the simplest level-1 case).
We may want to know if there is uncertainty about where gene flow came from,
that is, what groups are its sisters. For this, we can use functions currently
in PhyloNetworks, that map support onto a reference network:
[`PhyloNetworks.treeedges_support`](@extref) and
[`PhyloNetworks.hybridclades_support`](@extref).

Below, we use these functions to map the support of more features (like hybrid sisters) onto our level-1 consensus network as a reference network.
Note that, unlike the consensus functions here, the results of
`treeedges_support` and `hybridclades_support` depend on which hybrid edges
are *major* or *minor*, in the input networks.
The support for an edge to be a "tree" edge is the proportion of input networks
in which this edge is retained after deleting all the minor hybrid edges are
deleted, to get the major tree of each input network.

```@repl
BSe_majortree, majortree = treeedges_support(netsample, net);
BSn, BSe, BSc, BSgam, BSedgenum = hybridclades_support(netsample, net);
```

```@example
R"svg"(figname("consensus_7taxa_abc_hybsisters.svg"), width=9, height=3) # hide
R"layout"([1 2 3]); # hide
R"par"(mar=[0,0,0.5,0]); # hide
plot(net, edgelabel=BSe_majortree);
R"mtext"("support for edges to be in major tree", side=3, line=-1);
R"mtext"("support are in %, not proportions!", side=1, line=-1);
plot(net, edgelabeladj=[.5,-.2], nodelabeladj=[.5,.5],
    edgelabelcolor="orangered", nodelabelcolor="deepskyblue4",
    edgelabel=select(BSe,[:edge, :BS_hybrid_edge]),
    nodelabel=select(BSn,[:hybridnode,:BS_hybrid_samesisters])
);
R"mtext"("support for sister-hybrid relationships", side=3, line=-1);
R"mtext"("blue: hybrid + both sisters", side=1, line=-1);
plot(net, edgelabeladj=[.5,-.2],
    edgelabel=filter(row->row[:BS_hybrid]>0, BSn)[!,[:edge,:BS_hybrid]]);
R"mtext"("support for hybrid clades", side=3, line=-1);
R"mtext"("b1 of hybrid origin in 20% nets", side=1, line=-1);
R"dev.off()"; # hide
nothing # hide
```

![support for major-tree edges and sister-relationships, 5 input networks on 7 taxa](../assets/figures/consensus_7taxa_abc_hybsisters.svg)

From the left panel, we see that (b1,c1,c2) and (c1,c2) are in fact clades
in the major tree of 100% of our input networks.

See the section about PhyloNetworks's [Network support](@extref PhyloNetworks)
for more details.
