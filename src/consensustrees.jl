"""
    SplitFreq{N}

Frequency of one non-trivial taxon block, typically considered as describing an
unrooted bipartition (or split) of the `N` taxa into P=2 parts.

It is non-trivial if each part has at least 2 taxa, that is: the block is
not empty, not full, does not contain a single taxon, or all but a single taxon.

The taxon block is represented by an `N`-tuple of booleans, where entry `i`
says whether taxon number `i` is in or out of this block.
The canonical taxon block used to describe an unrooted bipartition is the block
that does not contain the last taxon, such that its last entry `N` is false.
"""
struct SplitFreq{N}
    """bipartition: block 1 does *not* contain the last taxon, block 2 does.
    The bipartition is described by 1 tuple of size N for membership in block 1:
    `split[i]` is `true` is taxon `i` is in block 1, `false` if it's in block 2.
    """
    split::NTuple{N,Bool}
    "frequency of the bipartition. mutable: use freq and freq! to get/set this value."
    freq::Base.RefValue{Int}
end
splitstring(obj::SplitFreq) = splitstring(obj.split)
splitstring_names(obj::SplitFreq, taxa::Vector) = splitstring_names(obj.split, taxa)

Base.show(io::IO, obj::SplitFreq{N}) where {N} = print(io,
    "SplitFreq on $N taxa, taxa in split cluster: " * splitstring(obj) *
    ", frequency: $(freq(obj))")

freq(obj::SplitFreq) = obj.freq[]
freq!(obj::SplitFreq, n) = obj.freq[] = n
incrementfreq!(obj::SplitFreq) = obj.freq[] += 1

"""
    split_fromHmatrix(M, i, N)

Tuple of booleans from row `i` of the hardwired-cluster matrix `M` on `N` taxa,
considering the phylogeny as unrooted: 0/1 values are switched if necessary,
to make sure that the last entry is `false`.

See [`cluster_fromHmatrix`](@ref) for the rooted version.
"""
function split_fromHmatrix(hwm, row::Int, N::Int)
    res = (hwm[row,N+1] == 1 ?
        ntuple(j -> !Bool(hwm[row, j+1]), N) :
        ntuple(j ->  Bool(hwm[row, j+1]), N)  )
    return res
end
"""
    cluster_fromHmatrix(M, i, N)

Tuple of booleans from row `i` of the hardwired-cluster matrix `M` on `N` taxa,
considering the phylogeny as *rooted*: 0/1 values as they are, *not* switched.

See [`split_fromHmatrix`](@ref) for the unrooted version.
"""
cluster_fromHmatrix(hwm, row::Int, N::Int) = ntuple(j ->  Bool(hwm[row, j+1]), N)

iscompatible(b1::SplitFreq{N}, b2::SplitFreq{N}) where N =
    treecompatible(b1.split, b2.split) # in utils.jl

function splitcomplement(splitvec::AbstractVector{NTuple{N,Bool}}) where N
    isoutgroup(i) = !any(t[i] for t in splitvec)
    return ntuple(isoutgroup, N)
end

function add_split!(bpvec::Vector{SplitFreq{N}}, split) where N
    i = findfirst(bp -> bp.split == split, bpvec)
    if isnothing(i)
        push!(bpvec, SplitFreq{N}(split,Ref(1)))
    else
        incrementfreq!(bpvec[i])
    end
end

"""
    consensustree(trees::AbstractVector{PN.HybridNetwork};
                  rooted=false,
                  proportion=0,
                  supportaslength=false)

Consensus tree summarizing the bipartitions (or clades) shared by more than
the required `proportion` of input `trees`.
An `ArgumentError` is thrown if one input network is not a tree, or the list of
input trees is empty, or if the input trees do not all have the same tip labels.
Input trees are not modified.

Output: consensus tree as an object of type `HybridNetwork`.

The bipartition (or clade) support values are stored
- as edge length with option `supportaslength=true` (*not* by default).
- in the field `.y` of each internal edge (as external edges correspond to
  trivial bipartitions, which must necessarily be in all sampled trees).
  This field is used by `writenewick` to write edge support, with its option
  `support=true`. However, this `.y` field is internal, so it can be modified
  by other functions and should not be relied upon.

By default, input trees are considered unrooted, and bipartitions are considered.
Use `rooted=true` to consider all input trees as rooted, in which case clades
(rather than bipartitions) are used to build the output rooted consensus tree.

By default, the greedy consensus consensus is calculated: the tree is built from
the bipartitions (or clades) with the highest support, until no more can be added.
The majority-rule tree can be obtained by using `proportion=0.5`: it is built
only from bipartitions (or clades) present in more than 50% of the input trees.

assumptions and **warnings**:
- Input trees are assumed to have their edges correctly directed.
  If unsure, run `directedges!.(trees)` prior.
- Input trees should not have degree-2 nodes other than the root
  (nodes with 1 only parent and 1 child).
  If unsure, run `removedegree2nodes!.(trees, true))` to keep their root even
  of degree 2 or `removedegree2nodes!.(trees, false))` unroot them also.

fixit: make a future version summarize edge lengths in
input trees and store their average in the consensus tree.

# example

```jldoctest
julia> nwk = ["((c,d),((a1,a2),b));", "(((a2,a1),b),c,d);", "(((a1,a2),c),d,b);"];

julia> treesample = readnewick.(nwk);

julia> con = consensustree(treesample); writenewick(con, round=true, support=true)
"(c,d,(b,(a1,a2)::1.0)::0.667);"

julia> con = consensustree(treesample; rooted=true); # greedy consensus

julia> writenewick(con, round=true, support=true)
"((d,c)::0.333,(b,(a1,a2)::1.0)::0.667);"

julia> [e.number => round(e.y, digits=3) for e in con.edge if !isexternal(e)] # edge number -> support
3-element Vector{Pair{Int64, Float64}}:
 6 => 0.333
 7 => 0.667
 8 => 1.0

julia> con = consensustree(treesample; rooted=true, proportion=0.5); # majority-rule

julia> writenewick(con, round=true, support=true)
"(c,d,(b,(a1,a2)::1.0)::0.667);"

julia> con = consensustree(treesample; proportion=0.75, supportaslength=true) |> writenewick
"(b,c,d,(a2,a1):1.0);"
```
"""
function consensustree(
    trees::AbstractVector{PN.HybridNetwork};
    rooted::Bool=false,
    proportion::Number=0,
    supportaslength::Bool=false
)
    isempty(trees) &&
        throw(ArgumentError("consensustree requires at least one network"))
    all(n.numhybrids==0 for n in trees) ||
        throw(ArgumentError("consensustree requires input trees (without reticulations)"))
    if length(trees) == 1
        net = deepcopy(trees[1])
        rooted || suppressroot!(net) # requires PN v1.3
        for e in net.edge
            if !isexternal(e)
                e.y = 1.0
                if supportaslength e.length = 1.0; end
            end
        end
        return net
    end
    taxa = sort!(tiplabels(trees[1]))
    ntaxa = length(taxa)
    splitcounts = SplitFreq{ntaxa}[]
    for net in trees
        length(net.leaf) == ntaxa ||
            throw(ArgumentError("input trees do not share the same taxon set"))
        count_bipartitions!(splitcounts, net, taxa, rooted) # will enforce same taxon sets
    end
    ntrees = length(trees)
    consensus_bipartitions!(splitcounts, proportion, ntrees)
    tree = tree_from_bipartitions(taxa, splitcounts, ntrees, supportaslength)
    for n in tree.node # clear internal non-hybrid node names
        if !n.leaf && !n.hybrid
            n.name = ""
        end
    end
    return tree
end

"""
    count_bipartitions!(counts, net, taxa, rooted)

Count bipartitions in `net`, and add them to `counts`.
If a new bipartition is found, a new key is added to `counts` with value 1.
If `net` has a bipartition already present as a key in `counts`, then the
corresponding value is incremented by 1.
`net` is not modified.

By default, input trees are considered unrooted, and bipartitions are considered.
Use `rooted=true` to consider all input trees as rooted, in which case clades
(rather than bipartitions) are used to build the output rooted consensus tree.

Only non-trivial (not all 0s or all 1s) splits contribute to the count totals. 

If the tip labels in `net` do not match those in `taxa` (as a set), then an
error will be thrown indirectly (via `PhyloNetworks.hardwiredclusters`).
"""
function count_bipartitions!(
    counts::Vector{SplitFreq{N}},
    net::PN.HybridNetwork,
    taxa::Vector{String},
    rooted::Bool,
) where N
    if !rooted
        if length(getroot(net).edge) < 3
            net = deepcopy(net) # re-binds the variable 'net'
            suppressroot!(net)
        end
    end
    getsplit = (rooted ? cluster_fromHmatrix : split_fromHmatrix)
    hw_matrix = hardwiredclusters(net, taxa)
    for row_idx in axes(hw_matrix, 1)
        split = getsplit(hw_matrix, row_idx, N)
        istrivialsplit(split) && continue
        add_split!(counts, split)
    end
    return counts
end

"""
    consensus_bipartitions!(splitcounts::Vector{SplitFreq{N}},
        proportion::Number, numtrees::Number)

Filter `splitcounts` to keep only the bipartitions whose frequencies are greater
than `proportion × numtrees`, or equal to `numtrees` when `proportion` is 1.
Bipartitions with frequency weight over 50% must be compatible with each other.
Bipartitions are retained one by one, from most to least frequent, so long as
they are compatible with bipartititions previously kept.

The result can be passed to [`tree_from_bipartitions`](@ref) to
construct the associated consensus tree topology.
`proportion = 1` corresponds to the strict consensus tree,
`proportion = 0.5` corresponds to the majority-rule consensus tree and
`proportion = 0` to a greedy consensus tree.

Output: `splitcounts` modified, with some entries filtered out, and sorted
by frequency (from least to most frequent) if `proportion<0.5`.

Assumption: all counts are positive.

Note that we use the term "bipartitions" here, which means splits of the taxon
set when ignoring the root. But this function can equally be used for clusters
(or clades) obtained as an edge's descendants when knowing the root position.

# Example
```jldoctest
julia> const PS = PhyloSummaries; # to use internals with less typing

julia> bp = [(true,false), (false,false), (true,true)]; freq=Ref.([3,1,4]);

julia> splitcounts = [PS.SplitFreq(x,y) for (x,y) in zip(bp, freq)]
3-element Vector{PhyloSummaries.SplitFreq{2}}:
 SplitFreq on 2 taxa, taxa in split cluster: 1, frequency: 3
 SplitFreq on 2 taxa, taxa in split cluster: , frequency: 1
 SplitFreq on 2 taxa, taxa in split cluster: 1,2, frequency: 4

julia> PS.consensus_bipartitions!(splitcounts, 0.5, 4)
2-element Vector{PhyloSummaries.SplitFreq{2}}:
 SplitFreq on 2 taxa, taxa in split cluster: 1, frequency: 3
 SplitFreq on 2 taxa, taxa in split cluster: 1,2, frequency: 4

julia> splitcounts = [PS.SplitFreq(x,y) for (x,y) in zip(bp, freq)]; # reset

julia> PS.consensus_bipartitions!(splitcounts, 0, 4) # sorted
3-element Vector{PhyloSummaries.SplitFreq{2}}:
 SplitFreq on 2 taxa, taxa in split cluster: , frequency: 1
 SplitFreq on 2 taxa, taxa in split cluster: 1, frequency: 3
 SplitFreq on 2 taxa, taxa in split cluster: 1,2, frequency: 4
```
"""
function consensus_bipartitions!(
    splitcounts::Vector{SplitFreq{N}},
    proportion::Number, 
    numtrees::Number,
) where N
    threshold2 = proportion * numtrees
    if proportion ≈ 1 # strict consensus
        filter!(v -> freq(v) ≈ numtrees, splitcounts)
    elseif threshold2 > 0 # 0 for greedy consensus: frequent case
        filter!(v -> freq(v) > threshold2, splitcounts)
    end
    sort!(splitcounts, by=freq, rev=false)
    if proportion >= 0.5
        return splitcounts
    end
    threshold1 = numtrees / 2 # all above must be compatible
    # next: traverse 'splitcounts' in reverse because we will delete items
    nS = length(splitcounts)
    for j_cs in nS:-1:1
        candidates = splitcounts[j_cs]
        freq(candidates) > threshold1 && continue
        iscompat = true
        for j_ks in (j_cs+1):length(splitcounts)
            if !iscompatible(candidates, splitcounts[j_ks])
                iscompat = false
                break
            end
        end
        iscompat || deleteat!(splitcounts, j_cs)
    end
    return splitcounts
end

"""
    tree_from_bipartitions(taxa::Vector{String},
        clustercounts::Vector{SplitFreq{N}},
        ntrees::Number,
        supportaslength::Bool)

Construct a tree on `N` taxa from a compatible set of clusters, as a
`PhyloNetworks.HybridNetwork` object.
Each cluster is represented as a tuple key `b`, and is given a weight:
its value `clusters[b]` divided by `ntrees`.
For each cluster, a node `n` and its parent edge `e` are added to the tree,
whose descendant taxa is the set `taxa[i]` for indices `i` such that b[i] is true.

The cluster's weight is stored in `e.y`.
With option `supportaslength=true`, the cluster weight is also in `e.length`.

Assumptions:
- `taxa` is of length `N`: same size as each cluster
- the input clusters are pairwise tree-compatible, which is the
  condition for them to be the clusters of a valid rooted tree.

Used by: [`consensustree`](@ref)
"""
function tree_from_bipartitions(
    taxa::Vector{String},
    clustercounts::Vector{SplitFreq{N}},
    ntrees::Number,
    supportaslength::Bool
) where N
    net = startree(taxa)
    ni = Ref(N+2) # leaves: 1:n, root n+1. next: n+2
    ei = Ref(N+1)
    for cf in clustercounts
        add_clusteredge!(net, ni, ei, cf.split, freq(cf)/ntrees, supportaslength)
    end
    return net
end

"""
    add_clusteredge!(tree, ni, ei, cluster, weight, supportaslength;
        nowarning=false)

Add a node (numbered `ni`) and its parent edge (numbered `ei`) in `tree`,
whose set of descendants is `cluster`, described by a 0/1 vector.
The node index `ni` is decremented, and the edge `ei` is incremented.

Output: newly created edge if successful (whose child node is the newly
created node), otherwise `nothing` with a warning unless `nowarning=true`.

The edge & node addition is unsuccessul if the cluster is empty, or is the full
taxon set, or if it is already present in `tree`. Adding it a second time
would result in a new node of degree 2, which is avoided.

Assumptions:
- `tree` is a tree (not checked)
- `cluster` is compatible with `tree` (considered rooted)
- leaf `j`, whose cluster membership is `cluster[j]`, is the node number `j` in
  the tree, that is, is it node `n = tree.node[i]` such that `n.number` is `j`.

Algorithm: the tree is traversed *before* being modified, to find:
- Q1: which node `lca` should be the parent of the new node, and
- Q2: which children of `lca` should become children of the new node.

To do so: `.booln2` and `.booln3` are used to store, for each node
* `.booln2`: Does 1+ descendant(s) of node ∈ clade ?
* `.booln3`: Is {node's descendants} ⊆ clade ?

Then, the answer to Q1 is the lowest node such that `.booln2 && !.booln3`,
and the answer to Q2 is all `lca`'s children with `.booln2 && .booln3`.
"""
function add_clusteredge!(
    net::PN.HybridNetwork,
    ni::Base.RefValue{Int},
    ei::Base.RefValue{Int},
    bv::SplitTuple,
    weight::Number,
    supportaslength::Bool;
    nowarning::Bool=false
)
    if all(bv) || all(.!bv)
        nowarning || @warn("will skip trivial clade: $bv")
        return nothing
    end
    # solve Q1: find lowest node with .booln2 && !.booln3
    lca, children_i = _lca_newcluster(net, bv)
    if length(children_i) == 1
        nowarning || @warn "clade already in tree (or incompatible?): will do nothing"
        return nothing
    end
    length(children_i) > 0 ||
        error("would not connect the new clade node to any children")
    # create a new node and new edge and
    # solve Q2: find all the lca's children with .booln2 && .booln3
    newedge = _resolveclade_belowlca(net, lca, children_i, ni, ei, weight, supportaslength)
    return newedge
end

#= (lca, ci) where `lca` solves Q1 and `ci` solves Q2:
- lca: LCA of cluster v, if it is not already in tree; parent of LCA otherwise.
  uses .booln2: does a node has ≥1 descendants in v?
  and  .booln3: are the nodes' descendants all in v?
  LCA: lowest node with .booln2 && !.booln3
- ci: indices in lca.edge of children of LCA whose descendants are all in v:
  children with .booln2 && .booln3, to solve Q2.
  These indices are sorted in reverse.
=#
@inline function _lca_newcluster(net, v)
    node2clade_intersection_initialize(net, v)   # initialize .booln{2,3}
    node2clade_intersection_update(getroot(net)) # post-order
    lca = getroot(net)
    (lca.booln2 && !lca.booln3) ||
        error("incorrect clade-intersection data at root, or trivial 1...1 clade")
    while true
        foundlca = true
        for ce in lca.edge
            cn = getchild(ce)
            cn === lca && continue
            if cn.booln2 && !cn.booln3
                lca = cn
                foundlca = false
                break
            end
        end
        foundlca && break # of while loop
    end
    children_i = _lca_newcluster_children(lca)
    return lca, children_i
end
# solve Q2
@inline function _lca_newcluster_children(lca)
    children_i = Int[]
    for i in length(lca.edge):-1:1
        cn = getchild(lca.edge[i])
        if cn !== lca && cn.booln2
            cn.booln3 || error("cn.booln3 should have been true")
            push!(children_i, i)
        end
    end
    return children_i
end

#= create a new node and new edge below lca, whose descendants is the cluster
used to compute the nodes' .booln2 and .booln3.
To be called after _lca_newcluster()
=#
@inline function _resolveclade_belowlca(net, lca,ci, ni,ei, wgt, supportaslength)
    newnode = PN.Node(ni[],false)
    newnode.name = "_$(ni[])"
    ni[] += 1
    # new edge: store clade support in .y, and as edge length if desired
    elen = (supportaslength ? wgt : -1.0)
    newe = PN.Edge(ei[],elen,false,wgt,-1.0,1.0, # z=-1, gamma=1
        [newnode,lca], true, # ischild1 is true: to agree with node ordering
        true,-1,true,true,false)
    ei[] += 1
    # Q2: loop over lca's children with .booln2 && .booln3
    for i in ci
        ce = lca.edge[i] # don't add new edge to lca.edge yet
        # cn = getchild(ce)
        # cn !== lca && cn.booln2 && cn.booln3 || error("oops")
        deleteat!(lca.edge,i) # disconnect lca-ce, connect newnode-ce
        push!(newnode.edge, ce)
        ce.node[2] = newnode # replaces lca, because ischil1 true and lca was parent
        # ischild1=true still agrees with newnode being the parent
    end
    # now we can modify lca.edge and net
    push!(lca.edge, newe)
    push!(newnode.edge, newe)
    PN.pushEdge!(net, newe)
    PN.pushNode!(net, newnode)
    return newe
end

# assumes: bv[j] corresponds to the node numbered j: j=net.node[i].number
function node2clade_intersection_initialize(net, bv)
    for node in net.node
        if node.leaf
            node.booln2 = node.booln3 = bv[node.number]
        else
            node.booln2 = false # OR of children
            node.booln3 = true  # AND of children
        end
    end
end
# assumes that
# 1. leaves have been initialized: true if in the clade, false otherwise
# 2. edges are correctly directed (correct ischild1 to use getchild)
function node2clade_intersection_update(pn::PN.Node)
    pn.leaf && return(nothing)
    for ce in pn.edge # loop over 'c'hild 'e'dges
        cn = getchild(ce)
        cn === pn && continue # skip parent edge
        node2clade_intersection_update(cn)
        if cn.booln2 && !pn.booln2 # OR from all children
            pn.booln2 = true
        end
        if !cn.booln3 && pn.booln3 # AND from all children
            pn.booln3 = false
        end
    end
    return nothing
end
