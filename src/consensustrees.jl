"""
    consensustree(trees::AbstractVector{PN.HybridNetwork};
                  rooted=false,
                  proportion=0)

Consensus tree summarizing the bipartitions (or clades) shared by more than
the required `proportion` of input `trees`.
An `ArgumentError` is thrown if one input network is not a tree, or the list of
input trees is empty, or if the input trees do not all have the same tip labels.
Input trees are not modified.

Output: consensus tree as an object of type `HybridNetwork`, with the
bipartition (or clade) support values stored as edge length, and also in
more "hidden" internal fields that should not be relied upon as they may change:
in each edge `.y` value, and in each node `.fvalue` for the clade descendant
from that node. Note that storing support values at nodes is *unsafe* for
unrooted bipartitions, as the bipartition associated with the node's descendants
given some rooting may become *one* of node's child subclade after re-rooting.

fixit:
- in docs/manual: include an example to plot the network with support shown
  for each edge
- use edge lengths in input trees...

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

# example

```jldoctest
julia> nwk = ["((c,d),((a1,a2),b));", "(((a2,a1),b),c,d);", "(((a1,a2),c),d,b);"];

julia> treesample = readnewick.(nwk);

julia> con = consensustree(treesample); writenewick(con, round=true)
"(c,d,((a2,a1):1.0,b):0.667);"

julia> con = consensustree(treesample; rooted=true); # greedy consensus

julia> writenewick(con, round=true)
(((a2,a1):1.0,b):0.667,(d,c):0.333);

julia> consensustree(treesample; rooted=true, proportion=0.5) # majority-rule
HybridNetwork, Semidirected Network
7 edges
8 nodes: 5 tips, 0 hybrid nodes, 3 internal tree nodes.
tip labels: a1, a2, b, c, ...
"(c,d,(b,(a1,a2):1.0):0.667);"

julia> consensustree(treesample; proportion=0.75) |> writenewick
"(b,c,d,(a2,a1):1.0);"
```
"""
function consensustree(
    trees::AbstractVector{PN.HybridNetwork};
    rooted::Bool=false,
    proportion::Number=0,
)
    isempty(trees) &&
        throw(ArgumentError("consensustree requires at least one network"))
    all(n.numhybrids==0 for n in trees) ||
        throw(ArgumentError("consensustree requires input trees (without reticulations)"))
    if length(trees) == 1
        net = deepcopy(trees[1])
        rooted || suppressroot!(net) # requires PN v1.3
        #= PN v1.3 is currently in branch 'master'.
        do this in some environment (but not in PhyloSummaries' main folder!)
        pkg> add PhyloNetworks#master
        these complications will go away after v1.3.0 is registered. =#
        # to preserve original lengths, do *not* store support as edge lengths
        # fixit: this is inconsistent with the case of 2+ input trees
        for e in net.edge
            if !isexternal(e)
                e.y = 1.0
                getchild(e).fvalue = 1.0
            end
        end
        return net
    end
    taxa = sort!(tiplabels(trees[1]))
    # taxa_set = Set(taxa)
    splitcounts = Dictionary{BitVector,Int}()
    for net in trees
        length(net.leaf) == length(taxa) ||
            throw(ArgumentError("input trees do not share the same taxon set"))
        #= not needed: will be checked by hardwiredclusters
        for leaf in net.leaf
            leaf.name in taxa_set ||
                throw(ArgumentError("taxon $(leaf.name) not in taxon list"))
        end
        =#
        count_bipartitions!(splitcounts, net, taxa, rooted)
    end
    ntrees = length(trees)
    consensus_bipartitions!(splitcounts, proportion, ntrees)
    return tree_from_bipartitions(taxa, splitcounts, ntrees)
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
    counts::Dictionary{BitVector,Int},
    net::PN.HybridNetwork,
    taxa::Vector{String},
    rooted::Bool,
) 
    if !rooted
        if length(getroot(net).edge) < 3
            net = deepcopy(net) # re-binds the variable 'net'
            suppressroot!(net)
        end
    end
    hw_matrix = hardwiredclusters(net, taxa)

    taxa_cols = 2:(length(taxa) + 1)
    for row_idx in axes(hw_matrix, 1)
        split = tuple_from_clustervector(view(hw_matrix, row_idx, taxa_cols), rooted)
        isnothing(split) && continue
        set!(counts, split, get(counts, split, 0) + 1)
    end
    return counts
end

"""
    tuple_from_clustervector(cluster01vector, rooted)

Tuple of booleans `t` with `t[i]` true / false if `taxa[i]` does / does not
belong in the hardwired cluster vector `cluster01vector` of 0/1 integers;
or `nothing` if the cluster is trivial (all 0s or all 1s).

If `rooted` is false, then clusters are considered as bipartitions and the last
taxon is used as an outgroup with a `false` entry.
For example, clusters `0011` and `1100` represent the same bipartition, and
both would return tuple `(true,true,false,false)`.
"""
function tuple_from_clustervector(cluster01vector::AbstractVector, rooted::Bool)
    n1 = sum(cluster01vector)
    n0 = length(cluster01vector)-n1
    if n1 <= 1 || n0 <= 1 # trivial clusters
        return nothing
    end
    if !rooted && cluster01vector[end] == 1
        return BitVector(x==0 for x in cluster01vector)
    end
    return BitVector(x==1 for x in cluster01vector)
end

"""
    consensus_bipartition!(splitcounts::Dictionary{BitVector,Int},
        proportion::Number, numtrees::Number)

Filter dictionary `splitcounts` to keep only the entries whose count (value in
the dictionary) greater than `proportion × numtrees`. Bipartitions with weight
over 50% must be compatible with each other. Bipartitions are added one by one,
from most to least frequent, so long as they are compatible with bipartititions
previously kept.

The result can be passed to [`tree_from_bipartitions`](@ref) to
construct the associated consensus tree topology.
`proportion = 0.5` corresponds to the majority-rule consensus tree and
`proportion = 0` to a greedy consensus tree.

Output: `splitcounts` modified, with some entries filtered out, and sorted
by frequency if `proportion<0.5`.

Assumption: all counts are positive.

# Example
```jldoctest
julia> using Dictionaries

julia> splitcounts = dictionary([[true,false]=>3, [false,false]=>1, [true,true]=>4])
3-element Dictionary{Vector{Bool}, Int64}:
 Bool[1, 0] │ 3
 Bool[0, 0] │ 1
 Bool[1, 1] │ 4

julia> consensus_bipartition!(splitcounts, 0.5, 4)
2-element Dictionary{Vector{Bool}, Int64}:
 Bool[1, 0] │ 3
 Bool[1, 1] │ 4

julia> splitcounts = dictionary([[true,false]=>3, [false,false]=>1, [true,true]=>4]);

julia> consensus_bipartition!(splitcounts, 0, 4)
3-element Dictionary{Vector{Bool}, Int64}:
 Bool[1, 1] │ 4
 Bool[1, 0] │ 3
 Bool[0, 0] │ 1
```
"""
function consensus_bipartitions!(
    splitcounts::Dictionary{BitVector,Int},
    proportion::Number, 
    numtrees::Number,
) 
    threshold1 = max(0.5, proportion) * numtrees # all above must be compatible
    threshold2 = proportion * numtrees  # applies if proportion < 0.5
    if proportion >= 0.5
        filter!(v -> v > threshold1, splitcounts)
        return(splitcounts)
    end
    if threshold2 > 0 # 0 for greedy consensus: frequent case
        filter!(v -> v > threshold2, splitcounts)
    end
    sort!(splitcounts, rev=false) # works for a Dictionary, but not a Dict
    nsplits = 0
    # next: traverse 'splitcounts' in reverse because we will delete items
    # use the "lazy" Iterators.reverse instead of Base.reverse, to avoid copying
    splits_mostfrequent = Iterators.reverse(keys(splitcounts)) # create this object once only
    for (candidate_bp, freq) in Iterators.reverse(pairs(splitcounts))
        if freq > threshold1
            nsplits += 1
            continue
        end
        # freq > threshold2 ensured by previous filtering
        iscompat = true
        for (i_bp, bp) in enumerate(splits_mostfrequent) # up-to-date because no copy
            i_bp > nsplits && break # only compare with previous splits
            if !treecompatible(candidate_bp, bp)
                iscompat = false
                break
            end
        end
        if iscompat # then keep candidate bipartition
            nsplits += 1
        else
            delete!(splitcounts, candidate_bp)
        end
    end
    return splitcounts
end

"""
    treecompatible(a::BitVector, b::BitVector)

true / false if two clusters `a` and `b` are / are not tree-compatible.

If A is the cluster of descendants of `a` (with `true` entries in `a`) and
if B is the cluster of descendant of `b`, these 2 clusters are tree-compatible
if there exists some tree that has both clusters.
This can be checked by the condition: A∩B is empty, or A⊆B, or B⊆A.
"""
function treecompatible(a::BitVector, b::BitVector)::Bool
    @assert length(a) == length(b)
    inter = a .& b
    if !any(inter)
        return true
    end
    if inter == a || inter == b
        return true
    end
    return false
end

"""
    tree_from_bipartitions(taxa::Vector{String},
        clusters::Dictionary{BitVector,<:Number},
        ntrees::Number)

Construct a consensus tree from a compatible set of cluster, as a
`PhyloNetworks.HybridNetwork` object.
Each cluster is represented as a `BitVector` key `b`, and is given a weight:
its value `clusters[b]` divided by `ntrees`.
For each cluster, a node is added to the consensus tree, whose descendant
taxa is the set `taxa[i]` for indices `i` such that b[i] is true.
The cluster's weight is stored in the node's `.fvalue` ---which is fragile,
if clusters are to be considered as unrooted bipartitions. The weight is also
stored in fields of the node's parent edge: `.length` and `.y`.

Assumption: the input clusters are pairwise tree-compatible, which is the
condition for them to be the clusters of a valid rooted tree.

Used by: [`consensustree`](@ref)
"""
function tree_from_bipartitions(
    taxa::Vector{String},
    bipartitions::Dictionary{BitVector,<:Number},
    ntrees::Number,
)
    n = length(taxa)
    net = PN.HybridNetwork()
    root = PN.Node(-2,false) # root has number -2
    # do *not* push the root to net.node yet: push leaves first, to make
    # taxa[i] be the label of net.node[i], so later we can use bv[leaf.number]
    # leaves: numbered 1:n
    for (i,t) in enumerate(taxa)
        edge = PN.Edge(i,-1.0) # ischild1 is true by default. length=-1 for NA
        leaf = PN.Node(i,true,false, # true: leaf
            -1.,[edge],false,false,false,false,false,false,-1,nothing,-1,-1,t)
        PN.pushNode!(net, leaf)
        edge.node = [leaf, root] # to match ischild1 is true
        PN.pushEdge!(net, edge)
        push!(root.edge, edge)
    end
    PN.pushNode!(net, root)
    net.rooti = n+1 # n leaves listed first, root listed next in net.node

    # internal nodes: numbered -3,-4 etc., as done by readnewick
    node_counter = -3
    edge_counter = n+1
    for (bv,weight) in pairs(bipartitions)
        all(bv) || all(.!bv) && @warn("will skip trivial clade: $bv")
        #= *before* modifying the network, traverse it to find
        Q1. which node 'lca' should be the parent of the new node, and
        Q2. which children of 'lca' should become children of the new node.
        To do so: use .booln2 and .booln3 to store:
           '1+ descendant of node ∈ clade?'
           '{node's descendants} ⊆ clade ?' =#
        node2clade_intersection_initialize(net, bv)
        node2clade_intersection_update(getroot(net)) # post-order
        # solve Q1: find lowest node with .booln2 && !.booln3
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
        # create a new node and new edge
        newnode = PN.Node(node_counter,false)
        weight /= ntrees
        newnode.fvalue = weight # store clade support in the clade's MRCA
        node_counter -= 1
        # new edge: store clade support as edge length and in .y
        newe = PN.Edge(edge_counter,weight,false,weight,0.0,1.0,
            [newnode,lca], true, # ischild1 is true: to agre with node ordering
            true,-1,true,true,false)
        edge_counter += 1
        # solve Q2: find all the lca's children with .booln2 && .booln3
        nchildren = 0
        nedges = length(lca.edge)
        for i in nedges:-1:1 # delete elements in lca.edge from the end
            ce = lca.edge[i] # don't add new edge to lca.edge yet
            cn = getchild(ce)

            cn === lca && continue
            if cn.booln2 && cn.booln3 # then cn should become a child of newnode
                nchildren += 1
                deleteat!(lca.edge,i) # disconnect lca-ce, connect newnode-ce
                push!(newnode.edge, ce)
                ce.node[2] = newnode # replaces lca, because ischil1 true and lca was parent
                # ischild1=true still agrees with newnode being the parent
            end
        end
        nchildren > 0 ||
            error("could not connect the new clade node to any children")
        # now we can modify lca.edge and net
        push!(lca.edge, newe)
        push!(newnode.edge, newe)
        PN.pushEdge!(net, newe)
        PN.pushNode!(net, newnode)
    end
    return net
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
