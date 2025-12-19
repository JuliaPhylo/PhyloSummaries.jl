#= treeofblobs
warn if one exit node has 2 edge
ig exit node does not have two try to match the cycle ordering
add a note saying that in the future we will add functionality to maintain circular order if the net is level 1
=#

#TODO: parse in circular order and store info that way. check funcs after canonicalizepartition

#= In a network with N leaves, the partition of taxa (leaves) defined by a blob
is represented as a tuple of N integers.
Taxa sharing the same integer are in the same taxon block (part of the partition).
This immutable type can serve as key is sets or dictionaries.
=#
const BlobSplit = NTuple{N,Bool} where N

"""
    BlobFreq{N,P}

Frequency of one blob partition, and frequency of its circular orders and
taxon blocks whose parent in the blob is a hybrid node.

- The partition associated with a blob in a network is the partition into
  taxon blocks from the connected component of the network after the blob
  (nodes and edges) is removed from the network.
- Different blobs with the same partition can have different circular orders
  of their taxon blocks. A blob does not necessarily admits a circular order.
- For each part in the partition, this part is a "hybrid" for this blob if it is
  adjacent to the blob at a hybrid node.

`P` is the number of parts (taxon blocks) in the partition.
`N` is the number of leaves (taxa) in the network.
"""
struct BlobFreq{N,P}
    "blob partition: P parts, each described as a tuple of of size N"
    partition::NTuple{P,NTuple{N,Bool}}
    "frequency of the blob partition. mutable: use freq and freq! to get/set this value."
    freq::Base.RefValue{Int}
    "frequencies of the different circular orders for blobs with this partition"
    circorder::Dict{NTuple{P,Int},Int}
    "frequencies of the different hybrid parts for blobs with this partition"
    hybrid::Dict{Int,Int}
    "bipartition (P=2) from a cut-edge, non-redundant with interesting blob partitions"
    cut::Bool
end
freq(obj::BlobFreq) = obj.freq[]
function freq!(obj::BlobFreq, n)
    obj.freq.x = n
    return obj
end
function incrementfreq!(obj::BlobFreq)
    obj.freq.x += 1
    return obj
end

"""
    consensus_treeofblobs(networks; proportion=0, minimumblobdegree=4)

Consensus tree summarizing the partitions of "interesting" blobs
(nodes of degree ≥ 4 in the blob tree) and the non-redundant bipartitions
(cut-edges connecting non-interesting blobs)
that are shared by more than the required `proportion` of input `trees`.
An error is thrown if the list of input networks is empty, or if the
input networks do not all have the same tip labels.

An "interesting" blob in an input network N is a non-trivial blob
(with at least one hybrid node) of degree m ≥ 4 by default.
The degree of a blob is the number of cut edges it is adjacent to,
and also the degree of the associated node in N's tree of blobs.
Setting `minimumblobdegree` to 3 or 2 will cause non-trivial blobs
to be considered "interesting" even if their corresponding node in N's
tree of blobs is of degree 3 or 2.

By default, a greedy consensus consensus is calculated.
The majority-rule tree can be obtained by using `proportion=0.5`,
and the strict consensus using `proportion=1`.

Assumptions:
- Each blob (2-edge connected component) is a biconnected component.
  This is true in binary phylogenies, or more generally when each tree node
  has at most 2 children and each hybrid node has a single child.
- fixit: polytomies...

todo:
- create another new function `consensus_level1network` similar to
  `consensus_treeofblobs`, that returns a network built with the same blobs as
  the consensus tree of blobs, but also resolves each blob with a level-1 cycle.
"""
function consensus_treeofblobs(
    networks::AbstractVector{PN.HybridNetwork};
    proportion::Number=0,
    minimumblobdegree::Int=4,
)
    isempty(networks) &&
        throw(ArgumentError("No input networks: cannot get a consensus"))
    minimumblobdegree ∈ (2,3,4) ||
        throw(ArgumentError("minimumblobdegree should be 2, 3 or 4, not $(minimumblobdegree)"))
    taxa = sort(tiplabels(networks[1]))
    blobarray = count_blobpartitions(networks, taxa, minimumblobdegree)
    # fixit: build consensus tree from blobarray using the proportion argument, and return it
end

"""
    count_blobpartitions!(networks, taxa, minBdegree)

Vector of [`BlobFreq{N}`](@ref) objects with `N` being the number of taxa,
with one entry for each blob partition: partition of all taxa into taxon blocks.
Each part, or taxon block, in each partition, is represented by a 0/1 tuple
with 1 at index `i` to indicate if `taxa[i]` is part of the block (0 if not).

Each object also has information about the frequency of the blob partition
in `networks`, the frequency of each circular order for that partition,
and frequency with which each taxon blob is hybrid for the blob.

A blob is a 2-edge connected component. In a non-binary network, a blob may
be the union of several biconnected components.

Special cases:
- Cut edges are trivial biconnected components, associated with a bipartition.
  A cut edge in some input network N contributes a `BlobFreq` entry if it is
  *not redundant* with a non-trivial blob of N, that is, if one of its two
  taxon blocks is not also a taxon block of a non-trivial blob in N.
- fixit: polytomies...

Note: the node field `.booln5` is used for intermediate calculations.

See also: [`consensus_treeofblobs`](@ref)
"""
function count_blobpartitions(
    networks::AbstractVector{PN.HybridNetwork},
    taxa::AbstractVector{<:String},
    minBdegree::Int,
)
    ntaxa = length(taxa)
    all(n.numtaxa == ntaxa for n in networks) ||
        throw(ArgumentError("input networks have different numbers of taxa"))
    # hardwiredclusters will error if different taxon sets
    blobarray = BlobFreq{ntaxa}[]
    for net in networks
        count_blobpartitions!(blobarray, net, taxa, minBdegree)
    end
    return blobarray
end

"""
    count_blobpartitions!(blobarray, net, taxa, minBdegree)
    count_blobpartitions!(blobarray, visitedbcc, net, taxaindex, minBdegree, blob, bidx, hwmatrix, edgemap)

Helper for [`count_blobpartitions`](@ref).
Update the entry of `blobarray` corresponding to the partition defined by
all blobs in one network `net`, or one `blob` in one network `net`.
If a new blob partition if found, that was absent from `blobarray`,
a new entry is created in the `blobarray` vector.

todo:
- double-check how we want to treat tree polytomies: 1 trivial blob with
  a single node, 0 hybrids, so corresponding to 0 biconnected component.
  Currently: "trivial" so "not interesting", so don't store their partition.
  Only store the bipartitions of adjacent non-redundant cut-edges.
- track which non-trivial blobs are "interesting"
- store non-redundant bipartitions: from trivial biconnected components not
  adjacent to an "interesting" blob already stored.
- use the minBdegree option: handle non-trivial blobs depending on their degree,
  and the downstream effect on which cut-edges are / are not redundant.
"""
function count_blobpartitions!(
    blobarray::Vector{BlobFreq{N}}, # shared number of taxa
    net::PN.HybridNetwork,
    taxa::AbstractVector{<:String},
    minBdegree::Int,
) where N
    taxaindex = Dict(t => i for (i, t) in pairs(taxa))
    PN.process_biconnectedcomponents!(net)
    hwmatrix = hardwiredclusters(net, taxa)
    edgemap = Dict(hwmatrix[i,1] => i for i in axes(hwmatrix,1))
    # add hybrid edges not listed in hwmatrix: because all partner edges on 1 row
    for h in net.hybrid
        pj = findfirst(e -> ischildof(h,e) && haskey(edgemap, e.number), h.edge)
        isnothing(pj) && error("no parent edge in hwcluster matrix for hybrid $(h.number)")
        pn = h.edge[pj].number
        for (j,e) in pairs(h.edge)
            (ischildof(h,e) && (j != pj)) || continue
            haskey(edgemap, e.number) &&
                error("2+ parent edges in hwcluster matrix for hybrid $(h.number)")
            edgemap[e.number] = pn
        end
    end
    visitedbcc = Set{Int}() # bicomponent already visited, from an earlier same blob
    for v in net.hybrid  # initialize before traversing all blobs
        v.booln5 = false # used by count_blobpartitions!
    end
    # pre-order traversal of biconnected components
    for (bidx, bc) in enumerate(net.partition)
        bidx ∈ visitedbcc && continue # bicomponent was already traversed
        PN.istrivial(bc) && continue
        # fixit (above): trivial blobs should not be ignored. Such a blob has
        # a single edge e=uv. It should be stored as a non-redundant bipartition
        # if both u and v are in blobs that are not "interesting" (not stored).
        count_blobpartitions!(blobarray, visitedbcc, net, taxaindex, minBdegree,
            bc, bidx, hwmatrix, edgemap)
    end
    return nothing
end
function count_blobpartitions!(
    blobarray::Vector{BlobFreq{N}}, # shared number of taxa
    visitedbcc::Set{Int}, # bicomponent reached & visited from an earlier one
    net::PN.HybridNetwork,
    taxaindex::Dict{String,Int},
    minBdegree::Int,
    blob::PN.Partition,
    bidx::Int,
    hwmatrix::AbstractMatrix,
    edgemap::Dict{<:Integer,<:Integer},
) where N
    blobdegree = Ref(0)
    splits, hybrids = blobtaxonsetpartition!(visitedbcc, blobdegree,
        net, blob, bidx, edgemap, hwmatrix, taxaindex, N)
    if blobdegree[] < minBdegree
        # fixit: store neighboring cut edges, because not redundant
        return blobarray
    end
    isempty(splits) && return blobarray
    # fixit: error if no split. Why ignore?
    isempty(hybrids) && return blobarray
    # fixit: error if no hybrid. Why ignore?
    partition = Tuple(splits)
    nparts = length(partition)
    # check if this partition already exists in blobarray
    # currently assuming only level 1
    matchidx, idxmap = findmatchingblob(blobarray, splits)
    if matchidx == -1
        # new blob
        defaultorder = ntuple(identity, nparts)
        circorder = Dict{typeof(defaultorder),Int}()
        circorder[defaultorder] = 1
        hybridmap = Dict{Int,Int}()
        hybridmap[hybrids[1]] = 1
        newblob = BlobFreq{N,nparts}(
            partition,
            Ref(1),
            circorder,
            hybridmap,
            false,
        )
        push!(blobarray, newblob)
    else
        # existing blob, increment freq
        bf = blobarray[matchidx]

        @assert length(hybrids) == 1 "expected a single hybrid index per blob"
        hybridpos = hybrids[1]

        # update hybrid frequency according to canonical partition slot
        canonhidx = idxmap[hybridpos]
        bf.hybrid[canonhidx] = get(bf.hybrid, canonhidx, 0) + 1

        # canonicalize circular orders using the split matched to partition entry 1
        startidx = findfirst(==(1), idxmap)
        startidx === nothing && return

        circorderkey, reversekey = canonicalorders(idxmap, startidx, hybridpos)
        isempty(circorderkey) && return
        if haskey(bf.circorder, circorderkey)
            bf.circorder[circorderkey] += 1
        elseif haskey(bf.circorder, reversekey)
            bf.circorder[reversekey] += 1
        else
            bf.circorder[circorderkey] = 1
        end
        incrementfreq!(bf)
    end
end

"""
    findmatchingblob(blobs, splits)

Return the index of the blob in `blobs` whose partition matches `splits`,
along with the permutation vector `idxmap`. Returns `(-1, Int[])` if no
match is found.
"""
function findmatchingblob(blobs::Vector{BlobFreq{N}}, splits::AbstractVector) where N
    for (i, blob) in pairs(blobs)
        partition = blob.partition
        length(partition) == length(splits) || continue

        idxmap = Vector{Int}(undef, length(splits))
        used = falses(length(partition))
        equalblob = true

        for (k, s) in pairs(splits)
            pos = findfirst(isequal(s), partition)
            if pos === nothing || used[pos]
                equalblob = false
                break
            end
            idxmap[k] = pos
            used[pos] = true
        end

        equalblob && return i, idxmap
    end
    return -1, Int[]
end

"""
    canonicalorders(idxmap, startidx, hybridpos)

Return the clockwise and counter-clockwise permutations of `idxmap`,
anchored at the split corresponding to the first partition entry, and
oriented relative to `hybridpos`.
"""
function canonicalorders(idxmap::AbstractVector{Int}, startidx::Int, hybridpos::Int)
    n = length(idxmap)
    if n == 0 || startidx < 1 || startidx > n ||
        hybridpos < 1 || hybridpos > n
        return (), ()
    end
    order = Int[]
    if startidx > hybridpos
        append!(order, idxmap[startidx:end])
        append!(order, reverse(idxmap[1:hybridpos]))
        if hybridpos + 1 <= startidx - 1
            append!(order, idxmap[(hybridpos+1):(startidx-1)])
        end
    else
        append!(order, idxmap[startidx:hybridpos])
        if hybridpos < n
            append!(order, reverse(idxmap[(hybridpos+1):n]))
        end
        if startidx > 1
            append!(order, idxmap[1:startidx-1])
        end
    end
    forward = Tuple(order)
    if length(order) > 1
        revseq = Vector{Int}(undef, length(order))
        revseq[1] = order[1]
        revseq[2:end] .= reverse(order[2:end])
        backward = Tuple(revseq)
    else
        backward = forward
    end
    return forward, backward
end


"""
    blobtaxonsetpartition!(visitedbcc, blobdegree, net, blob, bidx, edgemap, hwmatrix, taxaindex, ntaxa)

Depth-first traversal of a blob `B` that collect its taxon blocks, which form a
partition of the full set of `N` taxa. Each taxon block (or split) corresponds
to a cut-edge `uv` adjacent to the blob, with `u ∈ B` in the blob and `v ∉ B`.
The taxon block for `uv` is reprented by an `N`-tuple of 0/1 values with 1 at
index `i = taxaindex[label]` if the taxon named `label` is a descendant of
`uv`, and 0 otherwise.
If `u ∈ B` is a hybrid node, then the taxon block associated with `uv` is
considered "hybrid" for this blob.

Output: `(splits,hybrid)` where
- `splits` is the blob's partition as a tuple of N-tuples, and
- `hybrid` is the (sorted) vector of indices in `splits`, of taxon blocks
  that are hybrid for the blob.

Also:
- `visitedbcc` may be modified. If several biconnected components are part
  of the same blob starting at bicomponent indexed `bidx` (in `net.partition`),
  then the indices of these other bicomponents are added to `visitedbcc`.
- `blobdegree` is incremented by the number of taxon blocks found.

Warning: field `.booln5` is used to track if a hybrid node has (or not) been
visited yet. Should be initialized earlier.
"""
function blobtaxonsetpartition!(
    visitedbcc::Set{Int},
    blobdegree::Base.RefValue{Int},
    net::PN.HybridNetwork,
    bicomp::PhyloNetworks.Partition, # bicomponent: a blob may contain several
    bidx::Int, # bicomponent index: in net.partition, stored in e.inte1
    edgemap::Dict{<:Integer,<:Integer},
    hwmatrix::AbstractMatrix,
    taxaindex::Dict{String,Int},
    ntaxa::Int,
)
    entryidx  = PN.entrynode_preindex(bicomp)
    entrynode = net.vec_node[entryidx]
    splits = NTuple{ntaxa,Bool}[]
    hybrids = Int[]
    # traverse the blob starting at entrynode of the biconnected component
    blobtaxonsetpartition!(splits, hybrids, visitedbcc, blobdegree,
        entrynode, bidx, edgemap, hwmatrix, taxaindex, net)
    return splits, hybrids
end

"""
    blobtaxonsetpartition!(splits, hybrids, visitedbcc, blobdegree, entrynode, bidx, edgemap, hwmatrix, taxaindex, net)

Helper for `blobtaxonsetpartition` to accumulate entries in `splits` and `hybrids`.
The node field `.booln5` is used, to visit each hybrid node only once.
"""
function blobtaxonsetpartition!(
    splits::Vector{NTuple{N,Bool}},
    hybrids::Vector{Int},
    visitedbcc::Set{Int},
    blobdegree::Base.RefValue{Int},
    node::PN.Node,
    bidx::Int,
    edgemap::Dict{<:Integer,<:Integer},
    hwmatrix::AbstractMatrix,
    taxaindex::Dict{String,Int},
    net::PN.HybridNetwork,
) where N
    if node.hybrid # only visit a hybrid once
        node.booln5 && return nothing
        node.booln5 = true
    end
    taxacols = 2:(N + 1)
    # gather splits from 'node' *before* continuing the traversal
    ne_sp = 0 # number of edges giving a split: child cut-edges
    ne_go = 0 # number of edges to go through next: child edges in blob
    atotherBC = false
    bcc = net.partition
    for e in node.edge
        isparentof(node,e) || continue
        if e.inte1 == bidx # skip e if e ∈ this bicomponent
            ne_go += 1
            continue
        end
        if !PN.istrivial(bcc[e.inte1]) # e non-cut: e in ≠ non-trivial bicomp
            atotherBC = true           # can happen in non-binary net
            ne_go += 1
            continue
        end
        ne_sp += 1
        if getchild(e).leaf # trivial split, not in hwmatrix, but needed here
            i0 = taxaindex[getchild(e).name]
            split = Tuple(i == i0 for i in 1:N)
        else
            rowidx = get(edgemap, e.number, nothing)
            isnothing(rowidx) && error("unmapped non-external edge $(e.number)")
            # if good reason, then: split = descendants_bitvec(e, taxaindex)
            # split = Tuple(view(hwmatrix, rowidx, taxacols)) # tuple of Int
            split = Tuple(Bool(hwmatrix[rowidx,i]) for i in taxacols)
        end
        push!(splits, split)
        if node.hybrid
            push!(hybrids, length(splits))
        end
    end
    # fixit: add 1 more partition for the complement of all others (associated with entry node).
    # not needed if the root is the blob's entry node (= when bidx==1), needed otherwise.
    if ne_sp > 1 || atotherBC
        msg = "non-binary articulation node $(node.number): if its blob has a circular order, is is" *
            (ne_sp > 1 ? "" : " probably") * " not unique."
        @warn msg
    end
    blobdegree[] += ne_sp
    ne_go == 0 && return nothing # next loop not needed
    for e in node.edge
        isparentof(node,e) || continue
        ei = e.inte1
        if ei == bidx
            blobtaxonsetpartition!(splits, hybrids, visitedbcc, blobdegree,
                PN.getchild(e), bidx, edgemap, hwmatrix, taxaindex, net)
        elseif atotherBC && !PN.istrivial(bcc[ei])
            push!(visitedbcc, ei)
            otherBCentry = net.vec_node[PN.entrynode_preindex(bcc[ei])]
            otherBCentry.booln5 = false # reset: will go down (not back up)
            blobtaxonsetpartition!(splits, hybrids, visitedbcc, blobdegree,
                otherBCentry,   ei,   edgemap, hwmatrix, taxaindex, net)
        end
    end
    return nothing
end

# fixit: delete below? no longer used.
"""
    descendants_bitvec(edge::PN.Edge, taxaindex::Dict{String,Int})

Tuple `d` indicating which leaves are descendants of `edge`: `d[i]=true` if
a taxon "name" is a descendant of `edge`, `d[i]=false` otherwise,
with `i=taxaindex["name"]`.

Similar to `PhyloNetworks.descendants`, except that:
- the output contains `ntaxa` true/false values ordered according to `taxaindex`
  (`ntaxa` being the number of taxa)
  rather than the list of descendant node numbers, and
- only leaves are considered (`PhyloNetworks.descendants` can tell which
  internal nodes are descendants).
"""
function descendants_bitvec(edge::PN.Edge, taxaindex::Dict{String,Int})
    visited = Int[]
    splitvec = zeros(Bool, length(taxaindex))
    descendants_bitvec!(splitvec, visited, edge, taxaindex)
    return Tuple(splitvec)
end
function descendants_bitvec!(
    splitvec::Vector,
    visited::Vector{Int},
    edge::PN.Edge,
    taxaindex::Dict{String,Int}
)
    n = getchild(edge)
    if n.hybrid # only need to check previous visits for hybrid nodes
        n.number in visited && return nothing
        push!(visited, n.number)
    end
    if n.leaf
        splitvec[taxaindex[n.name]] = 1
    end
    for ce in n.edge
        if isparentof(n, ce)
            descendants_bitvec!(splitvec, visited, ce, taxaindex)
        end
    end
    return nothing
end
