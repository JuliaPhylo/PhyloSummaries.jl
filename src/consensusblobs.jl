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

Frequency of one non-trivial blob partition, and frequency of its circular orders
and taxon blocks whose parent in the blob is a hybrid node.

- The partition associated with a blob in a network is the partition into
  taxon blocks from the connected component of the network after the blob
  (nodes and edges) is removed from the network.
- Different blobs with the same partition can have different circular orders
  of their taxon blocks. A blob does not necessarily admits a circular order.
- For each part in the partition, this part is a "hybrid" for this blob if it is
  adjacent to the blob at a hybrid node.

`P` is the number of parts (taxon blocks) in the partition; 3 or more.
`N` is the number of leaves (taxa) in the network. 3 or more.
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
end

"""
    BipartFreq{N}

Same as [`BlobFreq`](@ref), but for a bipartition (P=2 parts) that is defined
by a cut-edge, not trivial, and not redundant with some interesting blob
(meaning: not adjacent to any interesting blob).
"""
struct BipartFreq{N}
    """bipartition: block 1 contains the last taxon, block 2 does not.
    The bipartition is described by 1 tuple of size N for membership in block 1:
    `split[i]` is `true` is taxon `i` is in block 1, `false` if it's in block 2.
    """
    split::NTuple{N,Bool}
    "frequency of the bipartition as a non-redundant. mutable: use freq and freq! to get/set this value."
    freq::Base.RefValue{Int}
end

freq(obj::Union{BipartFreq,BlobFreq}) = obj.freq[]
function freq!(obj::Union{BipartFreq,BlobFreq}, n)
    obj.freq[] = n
    return obj
end
function incrementfreq!(obj::Union{BipartFreq,BlobFreq})
    obj.freq[] += 1
    return obj
end

"""
    split_fromHmatrix(M, i, N)

Tuple of booleans from row `i` of the hardwired-cluster matrix `M` on `N` taxa,
considering the phylogeny as unrooted: 0/1 values are switched if necessary,
to make sure that the last entry is `false`.

Equivalent to `tuple_from_clustervector(M[i,2:(N+1)], false)`.
"""
function split_fromHmatrix(hwm, row::Int, N::Int)
    res = (hwm[row,N+1] == 1 ?
        ntuple(j -> !Bool(hwm[row, j+1]), N) :
        ntuple(j ->  Bool(hwm[row, j+1]), N)  )
    return res
end

"""
    consensus_treeofblobs(networks; proportion=0, minimumblobdegree=4)

Consensus tree summarizing the partitions of "interesting" blobs (nodes in
the tree of blobs) and the non-redundant bipartitions
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

Note that a node of degree 4 or more in the network's tree of blob may
correspond to a polytomy in N: a single node incident to m cut-edges, but
without any reticulation. These blobs are considered "non-interesting".
A cut-edge incident to such a polytomy is then non-redundant,
if the other blob it connects to is also non-intersting.

A chain of 2-blob leads to multiple cut-edge sharing the same bipartition.
This bipartition is counted only once (if non-trivial and non-redundant)
as if 2-blobs had been suppressed in the input network.

By default, a greedy consensus consensus is calculated.
The majority-rule tree can be obtained by using `proportion=0.5`,
and the strict consensus using `proportion=1`.

See also: [`count_blobpartitions!`](@ref)

todo:
- build the consensus tree from the vector of blobs and non-redundant bipartitions
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
    blobvec, bpvec = count_blobpartitions(networks, taxa, minimumblobdegree)
    # fixit: build consensus tree from these using the proportion argument, and return it
end

"""
    count_blobpartitions!(networks, taxa, minBdegree)

`(blob_vec, bipart_vec)` where `blob_vec` is a vector of [`BlobFreq{ntax}`](@ref)
object and `bipart_vec` is a vector of [`BipartFreq{ntax}`](@ref) objects,
`ntax` being the number of taxa.
All input networks must have the same set of `taxa`.

In `blob_vec`, each entry is for a non-trivial blob multi-partition:
a partition of all taxa into ≥ 3 taxon blocks.
Each part, or taxon block, in each partition, is represented by a 0/1 tuple
with 1 at index `i` to indicate if `taxa[i]` is part of the block (0 if not).
In `bipart_vec`, each entry is a bipartition (2 taxon blocs), counted if it
was not redundant with a non-trivial blob multi-partition.

Each object also has information about the frequency of the blob multi-partition
or non-redundant bipartition in `networks`. Each blob object also stores the
frequency of each circular order for that partition,
and the frequency with which each taxon blob is hybrid for the blob.

A blob is a 2-edge connected component. In a non-binary network, a blob may
be the union of several biconnected components.

A non-redundant bipartition comes from a cut-edge in an input network N.
A cut edge contributes and entry and is counted in a `bipart_vec` if it is
*not redundant* with a non-trivial blob of N, that is, if one of its two
taxon blocks is not also a taxon block of a non-trivial blob in N.

Side effects and internal fields:
- the node field `.intn1` stores 0 if the node is a singleton blob, and
  the node's blob index otherwise. This index is the index of the first
  bicomponent in the blob (which are pre-ordered).
- the node field `.booln5` is used for intermediate calculations.

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
    blobvec = BlobFreq{ntaxa}[]
    bpvec = BipartFreq{ntaxa}[]
    for net in networks
        count_blobpartitions!(blobvec, bpvec, net, taxa, minBdegree)
    end
    return blobvec, bpvec
end

"""
    count_blobpartitions!(blobs, biparts, net, taxa, minBdegree)

Helper for [`count_blobpartitions`](@ref).
Update the entries in the vector of `blobs` and in the vector of `biparts`
corresponding to the multi-partitions defined by all "interesting" blobs
and non-redundant bipartitions defined by cut-edges non-adjacent to some
interesting blob, in one network `net`.
If a new blob partition or non-redundant bipartition if found, that was absent
from `blobs` or `biparts` respectively, a new entry is created in the vector.

Notes:
- A "blob" here means a non-trivial blob with at least 1 hybrid node.
- A blob may be the union of 1 or more biconnected components
  (blobs partition nodes, bicomponents partition edges).
- A blob is "interesting" if its degree (number of adjacent cut edges)
  is at least `minBdegree`.
- Trivial biconnected components correspond to cut edges (and are not included
  in any blob, unlike non-trivial bicomponents). They are counted when they are
  *not redundant* with (not incident to) some "interesting" blob.
- A chain of 2-blob leads to multiple cut-edge sharing the same bipartition.
  This bipartition is counted only one (or not counted if it is trivial or
  adjacent to an interesting blob); as if 2-blobs had been suppressed.
"""
function count_blobpartitions!(
    blobvec::Vector{BlobFreq{N}}, # shared number of taxa
    bpvec::Vector{BipartFreq{N}},
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
    # initialize node fields used by count_blobpartitions! before traversal
    for v in net.hybrid
        v.booln5 = false # visited earlier during current bicomponent traversal?
    end
    for v in net.node
        v.intn1 = 0 # index of the node's blob, if non-trivial
    end
    # 'interesting' blobs will have blobdegree >= minBdegree
    blobdegree = zeros(Int, length(net.partition))
    # gather interesting blobs: pre-order traversal of biconnected components
    for (bidx, bc) in pairs(net.partition)
        bidx ∈ visitedbcc && continue # bicomponent was already traversed
        PN.istrivial(bc) && continue
        blobdegree[bidx] = count_blobpartitions!(blobvec, visitedbcc,
            net, taxaindex, minBdegree, bc, bidx, hwmatrix, edgemap)
    end
    # gather non-redundant cut-edges: from trivial bicomponents
    count_nonredundantbipartitions!(bpvec, blobdegree,
            net, taxaindex, minBdegree, hwmatrix, edgemap)
    return nothing
end

"""
    count_blobpartitions!(blobs, visitedbcc, net, taxaindex, minBdegree, blob, bidx, hwmatrix, edgemap)

Update the vector of `blobs` frequencies, and `visitedbcc` (to track
biconnected components already visited) for a single potentially interesting
non-trivial blob, starting from the `bidx`-th biconnected component in `net`.
"""
function count_blobpartitions!(
    blobvec::Vector{BlobFreq{N}}, # shared number of taxa
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
        return blobdegree[]
    end
    isempty(splits)  && error("non-trivial blob without any split.")
    isempty(hybrids) && error("non-trivial blob without any hybrid.")
    partition = Tuple(splits)
    nparts = length(partition)
    # check if this partition already exists in blobvec
    # currently assuming only level 1
    matchidx, idxmap = findmatchingblob(blobvec, splits)
    if isnothing(matchidx)
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
        )
        push!(blobvec, newblob)
    else
        # existing blob, increment freq
        bf = blobvec[matchidx]

        @assert length(hybrids) == 1 "expected a single hybrid index per blob"
        hybridpos = hybrids[1]

        # update hybrid frequency according to canonical partition slot
        canonhidx = idxmap[hybridpos]
        bf.hybrid[canonhidx] = get(bf.hybrid, canonhidx, 0) + 1

        # canonicalize circular orders using the split matched to partition entry 1
        startidx = findfirst(==(1), idxmap)
        startidx === nothing && return blobdegree[]

        circorderkey, reversekey = canonicalorders(idxmap, startidx, hybridpos)
        isempty(circorderkey) && return blobdegree[]
        if haskey(bf.circorder, circorderkey)
            bf.circorder[circorderkey] += 1
        elseif haskey(bf.circorder, reversekey)
            bf.circorder[reversekey] += 1
        else
            bf.circorder[circorderkey] = 1
        end
        incrementfreq!(bf)
    end
    return blobdegree[]
end

"""
    findmatchingblob(blobs::Vector{BlobFreq{N}}, splits) where N

Index in `blobs` and permutation to match the blob partition of `splits`.
Output:
- `(i,idxmap)` if `blobs[i]` matches `splits` using permutation vector `idxmap`,
  that is: `splits[k]` is `blobs[i].partition[idxmap[k]]` for all `k`.
- `(nothing, nothing)` if no match is found.

Assumption: the splits partition the full set of `N` taxa, and so do the
taxon blocks of each blob. In particular, splits are distinct from one another,
and taxon blocks from a given partition are also distinct.
"""
function findmatchingblob(
    blobs::Vector{BlobFreq{N}},
    splits::AbstractVector
) where N
    for (i, blob) in pairs(blobs)
        partition = blob.partition
        P = length(splits)
        length(partition) == P || continue
        idxmap = Vector{Int}(undef, P)
        # used = falses(length(partition)) ## removed because blocks are distincts
        equalblob = true
        for (k, s) in pairs(splits)
            pos = findfirst(isequal(s), partition)
            if pos === nothing
                equalblob = false
                break
            end
            idxmap[k] = pos
        end
        equalblob && return (i, idxmap)
    end
    return (nothing, nothing)
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

Warning: used internal node fields
- `.booln5` to track if a hybrid node has (or not) been visited yet.
- `.intn1` to track an index for the node's blob (which may contain more than 1
  bicomponent).
Both should be initialized earlier.
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
    entrynode.intn1 = bidx # blob index: index of first bicomponent
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
The node field `.booln5` is used to visit each hybrid node only once,
and `.intn1` is updated to store the blob index that visited nodes are in.
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
            split = ntuple(isequal(i0), N)
        else
            rowidx = get(edgemap, e.number, nothing)
            isnothing(rowidx) && error("unmapped non-external edge $(e.number)")
            # if good reason, then: split = descendants_bitvec(e, taxaindex)
            split = ntuple(i -> Bool(hwmatrix[rowidx,i+1]), N)
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
            cn = PN.getchild(e)
            cn.intn1 = node.intn1 # may be different from bidx
            blobtaxonsetpartition!(splits, hybrids, visitedbcc, blobdegree,
                cn, bidx, edgemap, hwmatrix, taxaindex, net)
        elseif atotherBC && !PN.istrivial(bcc[ei])
            push!(visitedbcc, ei)
            otherBCentry = net.vec_node[PN.entrynode_preindex(bcc[ei])]
            otherBCentry.booln5 = false # reset: will go down (not back up)
            otherBCentry.intn1 = node.intn1
            blobtaxonsetpartition!(splits, hybrids, visitedbcc, blobdegree,
                otherBCentry,   ei,   edgemap, hwmatrix, taxaindex, net)
        end
    end
    return nothing
end

"""
    count_nonredundantbipartitions!(bipart_vec, blobdegree, net, ...)

Gather non-redundant cut-edges from `net` and add them to `bipart_vec`
(or increment their frequency).
The field `.intn1` is used to know which blobs are adjacent to a cut edge,
and `blobdegree` to know if either of these blobs is "interesting"
"""
function count_nonredundantbipartitions!(
    bpvec::Vector{BipartFreq{N}},
    blobdegree::Vector{Int},
    net::PN.HybridNetwork,
    taxaindex::Dict{String,Int},
    minBdegree::Int,
    hwmatrix::AbstractMatrix,
    edgemap::Dict{<:Integer,<:Integer},
) where N
    inchain_store = Dict{Int,Bool}()
    inchain_leaf = Dict{Int,Union{String,Nothing}}()
    #= edge number => store?, for cut-edges at some end of a 2-blob chain:
    blobs of degree {d1,d2} = {2,d} with d≠2. If B is the blob of degree ≠2:
    false: B is interesting or d=1 (do not store the split)
    true:  B is not interesting (store the split if the other end agrees)
    =#
    for bc in net.partition
        PN.istrivial(bc)  || continue
        e = bc.edges[1] # edge from blob B1 -> blob B2
        n1 = getparent(e) # should be net.vec_node[p.cycle[1]]
        d1 = (n1.intn1 == 0 ? length(n1.edge) : blobdegree[n1.intn1])
        n2 = getchild(e)
        d2 = (n2.intn1 == 0 ? length(n2.edge) : blobdegree[n2.intn1])
        if d1 == 2 # edge in a 2-blob chain
            d2 == 2 && continue # inside the chain: skip
            inchain_store[e.number] = (1 < d2 < minBdegree) # bottom end: remember
            if d2 == 1
                inchain_leaf[e.number] = n2.name
            end
            continue # no storing decision yet
        end
        d2 == 1 && continue # trivial split (not end of 2-blob chain)
        if d2 == 2 # top end in a 2-blob chain
            inchain_store[e.number] = (1 < d1 < minBdegree) # top end: remember
            if d1 == 1 && n1.intn1 == 0
                inchain_leaf[e.number] = n1.name
            end
            continue # no storing decision yet
        end
        d1 == 1 && continue # can occur if root = leaf or ≠ LSA
        # by now, both d1>2 and d2>2
        (d1 < minBdegree && d2 < minBdegree) || continue # then e is non-redundant
        rowidx = get(edgemap, e.number, nothing)
        isnothing(rowidx) && error("unmapped non-external edge $(e.number)")
        split = ntuple(i-> Bool(hwmatrix[rowidx,i+1]), N)
        add_bipartition!(bpvec, split)
    end
    # add 0 or 1 biparts for each 2-blob chain: match the 2 end edges for each
    # 1. edges that have no entry in hwmatrix. trivial: don't store them
    while !isempty(inchain_leaf)
        e1, tax1 = pop!(inchain_leaf)
        ti = taxaindex[tax1]
        split = ntuple(isequal(ti), N)
        vsplitmatch(v) = v == split || all(v .!== split)
        row2 = findfirst(vsplitmatch, axes(hwmatrix,1))
        isnothing(row2) && error("blob of 2 chains without an internal edge")
        e2 = hwmatrix[row2, 1]
        haskey(inchain_store, e2) || error("2-blob chain: edge $e2 was not detected")
        pop!(inchain_store, e1)
        pop!(inchain_store, e2)
    end
    # 2. chains with internal edges at both ends: both in hwmatrix
    while !isempty(inchain_store)
        e1, store1 = pop!(inchain_store)
        haskey(edgemap, e1) || error("edge $e1 not in hw matrix")
        row1 = edgemap[e1]
        split = split_fromHmatrix(hwmatrix, row1, N)
        vsplitmatch(v) = v == split || all(v .!== split)
        isplitmatch(i) = i != row1 && splitmatch(view(hwmatrix,i,taxacols))
        row2 = findfirst(isplitmatch, axes(hwmatrix,1))
        isnothing(row2) && error("2-blob chains without 2 internal edges")
        e2 = hwmatrix[row2, 1]
        haskey(inchain_store, e2) || error("2-blob chain: edge $e2 was not detected")
        store2 = pop!(inchain_store, e2)
        if store1 && store2
            add_bipartition!(bpvec, split)
        end
    end
    return nothing
end

function add_bipartition!(bpvec::Vector{BipartFreq{N}}, split) where N
    i = findfirst(bp -> bp.split == split, bpvec)
    if isnothing(i)
        push!(bpvec, BipartFreq{N}(split,Ref(1)))
    else
        incrementfreq!(bpvec[i])
    end
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
