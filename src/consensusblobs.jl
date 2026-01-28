#TODO: parse in circular order and store info that way. check funcs after canonicalizepartition

"""
    BlobFreq{N,P}

Frequency of one non-trivial blob partition, and frequency of its circular orders
and taxon blocks whose parent in the blob is a hybrid node.

- The partition associated with a blob in a network is the partition into
  taxon blocks from the connected components of the network after the blob
  (nodes and edges) is removed from the network.
  Each part is represented by an `N`-tuple of booleans, where entry `i` says
  if taxon number `i` is in or out of this part.
- Different blobs with the same partition can have different circular orders
  of their taxon blocks. A blob does not necessarily admits a circular order.
- For each part in the partition, this part is a "hybrid" for this blob if it is
  adjacent to the blob at a hybrid node.

`P` is the number of parts (taxon blocks) in the partition; 3 or more.
`N` is the number of leaves (taxa) in the network. 3 or more.
"""
struct BlobFreq{N,P}
    """blob partition: P parts, each described as a tuple of of size N
    (an immutable type, so partitions can be keys in sets or dictionaries)"""
    partition::NTuple{P,NTuple{N,Bool}}
    "frequency of the blob partition. mutable: use freq and freq! to get/set this value."
    freq::Base.RefValue{Int}
    "frequencies of the different circular orders for blobs with this partition"
    circorder::Dict{NTuple{P,Int},Int}
    "frequencies of the different hybrid parts for blobs with this partition"
    hybrid::Dict{Int,Int}
end
partitionstring(obj::BlobFreq) = join(join.(findall.(obj.partition),","), "|")
Base.show(io::IO, obj::BlobFreq{N,P}) where {N,P} = print(io,
    "BlobFreq on $N taxa, $P blocks " * partitionstring(obj) *
    ", frequency $(freq(obj))" *
    ", $(length(obj.circorder)) circular orders, $(length(obj.hybrid)) hybrid blocks")
function Base.show(io::IO, ::MIME"text/plain", obj::BlobFreq{N,P}) where {N,P}
    println(io, "BlobFreq on $N taxa, partitioned into $P blocks")
    print(io, "taxon blocks: ")
    println(io, partitionstring(obj))
    print(io, "frequency: ")
    println(io, freq(obj))
    if !isempty(obj.circorder)
        print(io, "circular order of blocks => frequency: ")
        show(io, MIME"text/plain"(), obj.circorder)
    end
    print(io, "\nhybrid block => frequency: ")
    show(io, MIME"text/plain"(), obj.hybrid)
end


"""
    SplitFreq{N}

Frequency of one non-trivial taxon block, typically considered as describing an
unrooted bipartition (or split) of the `N` taxa into P=2 parts.

It is non-trivial if each part have at least 2 taxa, that is: the block is
not empty, not full, does not contain a single taxa, or all but a single taxon.

The taxon block is represented by an `N`-tuple of booleans, where entry `i`
says whether taxon number `i` is in or out of this block.
The canonical taxon block used to describe an unrooted bipartition is the block
that does not contain the last taxa, such that its last entry `N` is false.
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
splitstring(obj::SplitFreq) = join(findall(obj.split),",")
Base.show(io::IO, obj::SplitFreq{N}) where {N} = print(io,
    "SplitFreq on $N taxa, taxa in split cluster: " * splitstring(obj) *
    ", frequency: $(freq(obj))")

freq(obj::Union{SplitFreq,BlobFreq}) = obj.freq[]
function freq!(obj::Union{SplitFreq,BlobFreq}, n)
    obj.freq[] = n
    return obj
end
function incrementfreq!(obj::Union{SplitFreq,BlobFreq})
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

function splitcomplement(splitvec::AbstractVector{NTuple{N,Bool}}) where N
    isoutgroup(i) = !any(t[i] for t in splitvec)
    return ntuple(isoutgroup, N)
end

iscompatible(b1::SplitFreq{N}, b2::SplitFreq{N}) where N =
    treecompatible(b1.split, b2.split) # in utils.jl
iscompatible(b1::SplitFreq{N}, b2::BlobFreq{N}) where N =
    splitblobcompatible(b1.split, b2.partition)
iscompatible(b1::BlobFreq{N}, b2::BlobFreq{N}) where N =
    blobcompatible(b1.partition, b2.partition)
isredundantsplit(b1::SplitFreq{N}, b2::BlobFreq{N}) where N =
    isredundantsplit(b1.split, b2.partition)

"""
    consensus_treeofblobs(networks; proportion=0, minimumblobdegree=4)

Consensus tree summarizing the partitions of "interesting" blobs (nodes in
the tree of blobs) and the non-redundant bipartitions
(cut-edges connecting non-interesting blobs)
that are shared by more than the required `proportion` of input `trees`.
An error is thrown if the list of input networks is empty, or if the
input networks do not all have the same tip labels.

!!! note "This is an unrooted tree"
    This tree is to be considered as unrooted. It is built arbitrarily
    rooted at the last alphabetical taxon.
    Use `rootatnode!` or `rootonedge!` to re-root this tree
    given external knowledge of the outgroup (taxon or clade).

The support for a blob is the proportion of input networks that have a blob
with this partition. This is stored in the corresponding node's `.fvalue`.
The support for a bipartition as non-redundant is the proportion of input
networks that have this bipartition *not* adjacent to any interesting blob.
This is stored in the corresponding edge's field `.y`.
With option `supportaslength=true`, this is also stored in the edge's
`.length`. This option is *not* recommended and may be removed.

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
"""
function consensus_treeofblobs(
    networks::AbstractVector{PN.HybridNetwork};
    proportion::Number=0,
    minimumblobdegree::Int=4,
    supportaslength::Bool=false,
)
    isempty(networks) &&
        throw(ArgumentError("No input networks: cannot get a consensus"))
    minimumblobdegree ∈ (2,3,4) ||
        throw(ArgumentError("minimumblobdegree should be 2, 3 or 4, not $(minimumblobdegree)"))
    taxa = sort(tiplabels(networks[1]))
    blobvec, bpvec = count_blobpartitions(networks, taxa, minimumblobdegree)
    nnets = length(networks)
    filter_sort_compatible_partitions!(blobvec, bpvec, nnets, proportion)
    tob, _ = tree_from_blobpartitions(taxa, blobvec, bpvec, nnets, supportaslength)
    return tob
end

"""
    consensus_level1network(networks; proportion=0, minimumblobdegree=4)

fixit
"""
function consensus_level1network(
    networks::AbstractVector{PN.HybridNetwork};
    proportion::Number=0,
    minimumblobdegree::Int=4,
    supportaslength::Bool=false,
)
    isempty(networks) &&
        throw(ArgumentError("No input networks: cannot get a consensus"))
    minimumblobdegree ∈ (3,4) ||
        throw(ArgumentError("minimumblobdegree should be 3 or 4, not $(minimumblobdegree)"))
    taxa = sort(tiplabels(networks[1]))
    blobvec, bpvec = count_blobpartitions(networks, taxa, minimumblobdegree)
    hybdict = count_hybridclusters(blobvec) # before filtering blobs out
    nnets = length(networks)
    filter_sort_compatible_partitions!(blobvec, bpvec, nnets, proportion)
    net, fixit = tree_from_blobpartitions(taxa, blobvec, bpvec, nnets, supportaslength)
    edge2hyb = _hybridsupport_cladesintree(hybdict, tob, taxa, nnets)
    # fixit: expand each blob node into a cycle: most frequent circular order, most frequent compatible hybrid
    return net
end

"""
    count_blobpartitions!(networks, taxa, minBdegree)

`(blob_vec, bipart_vec)` where `blob_vec` is a vector of [`BlobFreq{ntax}`](@ref)
object and `bipart_vec` is a vector of [`SplitFreq{ntax}`](@ref) objects,
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
- `.inte1` set by `process_biconnectedcomponents!` is used (and not modified)
- `.intn1` stores 0 if a node is a singleton blob, and the node's blob index
  otherwise: index of the first bicomponent in the blob (which are pre-ordered).
- `.boole1` of edges, to visit hybrid nodes once and in "half" circular order.

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
    bpvec = SplitFreq{ntaxa}[]  # bipartitions, frequency: if non-redundant
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
    bpvec::Vector{SplitFreq{N}},
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
    # initialize fields used by count_blobpartitions! before traversal
    for e in net.edge
        e.boole1 = true # may be traversed?
    end
    for v in net.node
        v.intn1 = 0 # index of the node's blob. will remain 0 if {v} = trivial blob
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
    blobdegree = (bidx==1 ? Ref(0) : Ref(1)) # parent edge: 0 if root, 1 ow
    splits, hybrids, islevel1 = blobtaxonsetpartition!(visitedbcc, blobdegree,
        net, blob, bidx, edgemap, hwmatrix, taxaindex, N)
    if blobdegree[] < minBdegree
        return blobdegree[]
    end
    isempty(splits)  && error("non-trivial blob without any split.")
    isempty(hybrids) && error("non-trivial blob without any hybrid.")
    nparts = length(splits)
    # check if this partition already exists in blobvec
    matchidx, idxmap = findmatchingblob(blobvec, splits)
    if isnothing(matchidx) # add new blob to blobvec
        if islevel1
            defaultorder = ntuple(identity, nparts)
            #= fixit: 'identity' wrong if 2+ blocks on both sides of the hybrid.
            ihyb = hybrids[1]
            side1 = 1:ihyh; side2 = nparts:-1:(ihyb+1)
            defaultorder = Tuple(vcat(side1, side2))
            or something more efficient?
            downstream problems with changing the default order?
            Alternative: build the Tuple of splits with this order instead?
            co = vcat(side1, side2) # or collect(Iterators.flatten((side1, side2)))
            partition = ntuple(i -> split[co[i]], nparts)
            newblob = BlobFreq{N,nparts}(partition, Ref(1), cofreq, hybmap)
            =#
            cofreq = Dict(defaultorder => 1)
        else
            # level > 1: do not store circular order
            cofreq = Dict{NTuple{nparts,Int},Int}()
        end
        hybmap = Dict(hybrids[1] => 1)
        newblob = BlobFreq{N,nparts}(Tuple(splits), Ref(1), cofreq, hybmap)
        push!(blobvec, newblob)
    else # existing blob: increment frequencies of canonical partition slots
        bf = blobvec[matchidx]
        length(hybrids) == 1 || @warn "expected a single exit hybrid per blob"
        hybridpos = hybrids[1] # use the first to find circular order
        for k in hybrids       # but count all
            canonk = idxmap[k]
            bf.hybrid[canonk] = get(bf.hybrid, canonk, 0) + 1
        end
        # only calculate and store circular order for level-1 blobs
        if islevel1
            # canonicalize circular orders using the split matched to partition entry 1
            startidx = findfirst(==(1), idxmap)
            startidx === nothing && error("blob partition: 1 not found in idxmap. Expected Set(idxmap) == Set(1:nparts)")
            #fixit:
            # we should always have that Set() == Set(1:nparts)
            # also, could this help: indexin(1:nparts, idxmap)
            # first value = startidx
            circorderkey, reversekey = canonicalorders(idxmap, startidx, hybridpos)
            isempty(circorderkey) && error("blob partition: empty circular order key")
            if haskey(bf.circorder, circorderkey)
                bf.circorder[circorderkey] += 1
            elseif haskey(bf.circorder, reversekey)
                bf.circorder[reversekey] += 1
            else
                bf.circorder[circorderkey] = 1
            end
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
    elseif startidx == hybridpos
        append!(order, idxmap[startidx:end])
        append!(order, reverse(idxmap[1:hybridpos-1]))

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
If `u ∈ B` is a hybrid node incident to some exit edge `uv`, then the taxon
block associated with this edge is considered "hybrid" for this blob.

Output: `(splits,hybrid)` where
- `splits` is the blob's partition as a tuple of N-tuples, listed in a "half"
  circular order if the blob is level-1: from highest to lowest along one side
  then along the other side.
- `hybrid` is the vector of indices in `splits`, of taxon blocks
  that are hybrid for the blob.

Also:
- `visitedbcc` may be modified. If several biconnected components are part
  of the same blob starting at bicomponent indexed `bidx` (in `net.partition`),
  then the indices of these other bicomponents are added to `visitedbcc`.
- `blobdegree` is incremented by the number of taxon blocks found.

Warning: used internal fields
- `.inte1` set by `process_biconnectedcomponents!` is used (not modified)
- `.intn1` is modified to track to store the node's blob, if non-trivial
  (a blob may contain more than 1 bicomponent). This field should be
  initialized earlier (to a value ≤ 0).
- `.boole1` is modified to track if a hybrid edge may be traversed, still;
  should be initialized earlier (to true).
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
    # bloblevel: true if all biconnected components in this blob are level-1
    bloblevel = Ref(PN.getlevel(bicomp) <= 1)
    # traverse the blob starting at entrynode of the biconnected component
    blobtaxonsetpartition!(splits, hybrids, visitedbcc, blobdegree, bloblevel,
        entrynode, bidx, edgemap, hwmatrix, taxaindex, net)
    if bidx > 1 # entry ≠ root: add split from the unique entry cut-edge
        outgroupsplit = splitcomplement(splits)
        if any(outgroupsplit) # empty if entry is at or above LSA (≠ root)
            # first to get circular order later: treat a node before its children
            pushfirst!(splits, outgroupsplit)
            hybrids .+= 1
        end
    end
    return splits, hybrids, bloblevel[]
end

"""
    blobtaxonsetpartition!(splits, hybrids, visitedbcc, blobdegree, bloblevel, entrynode, bidx, edgemap, hwmatrix, taxaindex, net)

Helper for `blobtaxonsetpartition` to accumulate entries in `splits` and `hybrids`.
Internal fields:
- `.inte1` used but not modified: should store the index of an edge's
  biconnected component.
- `.intn1` updated to store the blob index that a visited node is in:
  index of the first biconnected component in this blob.
-`.boole1` set to `false` for all partners of an edge that will be traverered.

In `splits`, only the *children* taxon blocks are gathered, from exit cut-edges,
but recursively across all bicomponents in the blob (which may occur in
non-binary networks).
"""
function blobtaxonsetpartition!(
    splits::Vector{NTuple{N,Bool}},
    hybrids::Vector{Int},
    visitedbcc::Set{Int},
    blobdegree::Base.RefValue{Int},
    bloblevel::Base.RefValue{Bool},
    node::PN.Node,
    bidx::Int,
    edgemap::Dict{<:Integer,<:Integer},
    hwmatrix::AbstractMatrix,
    taxaindex::Dict{String,Int},
    net::PN.HybridNetwork,
) where N
    # gather splits from 'node' *before* continuing the traversal
    ne_sp = 0 # number of edges giving a split: child cut-edges
    ne_go = 0 # number of edges to go through next: child edges in blob
    atotherBC = false
    bcc = net.partition
    for e in node.edge
        isparentof(node,e) || continue
        e.boole1 || continue # hybrid node already visited
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
            split = ntuple(i -> Bool(hwmatrix[rowidx,i+1]), N)
        end
        push!(splits, split)
        if node.hybrid
            push!(hybrids, length(splits))
        end
    end
    if ne_sp > 1 || atotherBC
        msg = "non-binary articulation node $(node.number): if its blob has a circular order, is is" *
            (ne_sp > 1 ? "" : " probably") * " not unique."
        @warn msg
    end
    blobdegree[] += ne_sp
    ne_go == 0 && return nothing # next loop not needed
    for e in node.edge
        isparentof(node,e) || continue
        e.boole1 || continue
        ei = e.inte1
        if ei == bidx
            cn = PN.getchild(e)
            if e.hybrid # then partner edges should be skipped later
                for pe in cn.edge
                    if pe.hybrid && pe !== e && ischildof(cn,pe)
                        pe.boole1 = false
                    end
                end
            end
            cn.intn1 = node.intn1 # may be different from bidx
            blobtaxonsetpartition!(splits, hybrids, visitedbcc, blobdegree, bloblevel,
                cn, bidx, edgemap, hwmatrix, taxaindex, net)
        elseif atotherBC && !PN.istrivial(bcc[ei]) && !(ei ∈ visitedbcc)
            push!(visitedbcc, ei)
            otherBCentry = net.vec_node[PN.entrynode_preindex(bcc[ei])]
            otherBCentry.intn1 = node.intn1
            # check if other BCC is level > 1: if so, blob is not level-1
            if PN.getlevel(bcc[ei]) > 1
                bloblevel[] = false
            end
            blobtaxonsetpartition!(splits, hybrids, visitedbcc, blobdegree, bloblevel,
                otherBCentry,   ei,   edgemap, hwmatrix, taxaindex, net)
        end
    end
    return nothing
end
#Tstore the edge that connects blobs we dont count
"""
    count_nonredundantbipartitions!(bipart_vec, blobdegree, net, ...)

Gather non-redundant cut-edges from `net` and add them to `bipart_vec`
(or increment their frequency).
The field `.intn1` is used to know which blobs are adjacent to a cut edge,
and `blobdegree` to know if either of these blobs is "interesting"
"""
function count_nonredundantbipartitions!(
    bpvec::Vector{SplitFreq{N}},
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
        s1 = n1.intn1 == 0 || d1 < minBdegree # B1 not interesting: store e (perhaps)
        n2 = getchild(e)
        d2 = (n2.intn1 == 0 ? length(n2.edge) : blobdegree[n2.intn1])
        s2 = n2.intn1 == 0 || d2 < minBdegree
        if d1 == 2 # edge in a 2-blob chain
            d2 == 2 && continue # inside the chain: skip
            inchain_store[e.number] = s2 && d2 > 1 # bottom end: remember
            if d2 == 1
                inchain_leaf[e.number] = n2.name
            end
            continue # no storing decision yet
        end
        d2 == 1 && continue # trivial split (not end of 2-blob chain)
        if d2 == 2 # top end in a 2-blob chain
            inchain_store[e.number] = s1 && d1 > 1 # top end: remember
            if d1 == 1 && n1.intn1 == 0
                inchain_leaf[e.number] = n1.name
            end
            continue # no storing decision yet
        end
        d1 == 1 && continue # can occur if root = leaf or ≠ LSA
        # by now, both d1>2 and d2>2
        (s1 && s2) || continue # skip if redundant with B1 or B2
        row = get(edgemap, e.number, nothing)
        isnothing(row) && error("unmapped non-external edge $(e.number)")
        split = split_fromHmatrix(hwmatrix, row, N)
        add_bipartition!(bpvec, split)
    end
    # add 0 or 1 biparts for each 2-blob chain: match the 2 end edges for each
    # 1. edges that have no entry in hwmatrix. trivial: don't store them
    while !isempty(inchain_leaf)
        e1, tax1 = pop!(inchain_leaf)
        ti = taxaindex[tax1]
        split = ntuple(isequal(ti), N)
        vsplitmatch(v) = v == split || all(v .!== split)
        rows = findall(vsplitmatch, axes(hwmatrix,1))
        k = findfirst(r -> hwmatrix[r,1] != e1 && haskey(inchain_store, hwmatrix[r,1]), rows)
        isnothing(k) && error("blob of 2 chains without an internal edge")
        e2 = hwmatrix[rows[k], 1]
        haskey(inchain_store, e2) || error("2-blob chain: edge $e2 was not detected")
        pop!(inchain_store, e1)
        pop!(inchain_store, e2)
    end
    # 2. chains with internal edges at both ends: both in hwmatrix
    while !isempty(inchain_store)
        e1, store1 = pop!(inchain_store)
        row1 = get(edgemap, e1, nothing)
        isnothing(row1) && error("edge $e1 not in hw matrix")
        split = split_fromHmatrix(hwmatrix, row1, N)
        vsplitmatch(v) = v == split || all(v .!== split)
        isplitmatch(i) = i != row1 && vsplitmatch(view(hwmatrix,i,2:(N+1)))
        rows = findall(isplitmatch, axes(hwmatrix,1))
        k = findfirst(r -> haskey(inchain_store, hwmatrix[r,1]), rows)
        isnothing(k) && error("2-blob chains without 2 internal edges")
        e2 = hwmatrix[rows[k], 1]
        haskey(inchain_store, e2) || error("2-blob chain: edge $e2 was not detected")
        store2 = pop!(inchain_store, e2)
        if store1 && store2
            add_bipartition!(bpvec, split)
        end
    end
    return nothing
end

function add_bipartition!(bpvec::Vector{SplitFreq{N}}, split) where N
    i = findfirst(bp -> bp.split == split, bpvec)
    if isnothing(i)
        push!(bpvec, SplitFreq{N}(split,Ref(1)))
    else
        incrementfreq!(bpvec[i])
    end
end


"""
    filter_sort_compatible_partitions!(blobpartitions, bipartitions, nnets, proportion)

Filter out blob- and bi-partitions with frequency ≤ proportion;
sort each vector of `blobpartitions` and `bipartitions` by frequency
(from smallest to largest);
filter out any partition not compatible with another of higher frequency
(giving preference to a blob over a bipartition if tied).

Blob-compatibility is used, which reduces to tree-compatibility when both
partitions are bipartitions.

Bipartitions redundant with a retained blob are *not* filtered out.
"""
function filter_sort_compatible_partitions!(
    blobparts::Vector{BlobFreq{N}},
    biparts::Vector{SplitFreq{N}},
    nnets::Number,
    proportion::Number,
) where N
    threshold2 = proportion * nnets
    if proportion ≈ 1 # strict consensu
        filter!(v -> freq(v) ≈ nnets, blobparts)
        filter!(v -> freq(v) ≈ nnets, biparts)
    elseif threshold2 > 0 # 0 for greedy consensus: frequent case
        filter!(v -> freq(v) > threshold2, blobparts)
        filter!(v -> freq(v) > threshold2, biparts)
    end
    sort!(blobparts, by=freq, rev=false)
    sort!(biparts,   by=freq, rev=false)
    if proportion >= 0.5 # fixit: and if binary networks and if minBdegree >=4
        # then: no need to check for compatibilitiy
        return nothing
    end
    threshold1 = max(0.5, proportion) * nnets
    # filter for within-list compatibility
    for bparts in (blobparts, biparts)
        nB = length(bparts)
        for j_cb in nB:-1:1
            candidateb = bparts[j_cb]
            freq(candidateb) > threshold1 && continue
            iscompat = true
            for j_kb in (j_cb+1):length(bparts)
                if !iscompatible(candidateb, bparts[j_kb])
                    iscompat = false
                    break
                end
            end
            iscompat || deleteat!(bparts, j_cb)
        end
    end
    # filter for between-list compatibility
    nbb_j = length(blobparts) # Next BloB / BIpartition index to check:
    nbi_j = length(biparts)   # from most to least frequent
    nbb_f = (nbb_j>0 ? freq(blobparts[nbb_j]) : 0)
    nbi_f = (nbi_j>0 ? freq(biparts[nbi_j])   : 0)
    while nbb_j>0 || nbi_j>0
        # if equally frequent: favor keeping the blob
        if nbi_f > nbb_f # decide to keep or filter out the next bipart
            cb = biparts[nbi_j] # Candidate Bipartition
            keep = true
            for j_kb in (nbb_j+1):length(blobparts)
                kb = blobparts[j_kb] # Kept Blob
                if !iscompatible(cb, kb)
                    keep = false
                    break
                end
            end
            keep || deleteat!(biparts, nbi_j)
            nbi_j -= 1
            nbi_f = (nbi_j>0 ? freq(biparts[nbi_j]) : 0)
        else # decide to keep or filter out the next blob
            cb = blobparts[nbb_j]
            iscompat = true
            for j_kb in (nbi_j+1):length(biparts)
                if !iscompatible(biparts[j_kb], cb) # okay if redundant
                    iscompat = false
                    break
                end
            end
            iscompat || deleteat!(blobparts, nbb_j)
            nbb_j -= 1
            nbb_f = (nbb_j>0 ? freq(blobparts[nbb_j]) : 0)
        end
    end
    return nothing
end

"""
    tree_from_blobpartitions

Tree summarizing the input partitions: with a node for each input blob
and an edge (if not redundant) for each input bipartition.
This tree should be considered unrooted, in part because blob partitions
do not have root information (all their taxon blocks are listed).

A blob's weight is stored in the corresponding node's `.fvalue`.
A bipartition's weight is stored in the corresponding edge's field `.y`.
With option `supportaslength=true`, this weight is also stored in
the bipartition edge's `.length`.

Assumptions (not checked):
1. each bipartition is represented by the taxon block that does *not* contain
   the last taxon
2. bipartititions are not trivial: they separate at least 2 taxa from at 2 other
3. blobs are correct partitions, with P ≥ 3 parts at least
4. blob partitions and bipartitions are all *compatible* with each other:
   any 2 blobs from `blobpartitions` are blob-compatible;
   any 2 splits from `bipartitions` are tree-compatible; and
   any blob and split are blob-compatible.

Calls [`add_blobnode!`](@ref) and [`add_clusteredge!`](@ref), which use
internal fields `.booln2` and `.booln3`.
"""
function tree_from_blobpartitions(
    taxa::Vector{String},
    blobparts::Vector{BlobFreq{N}},
    biparts::Vector{SplitFreq{N}},
    nnets::Number,
    supportaslength::Bool
) where N
    N == length(taxa) || error("N ($N) != number of taxa $(length(taxa))")
    # fixit: initialize vector of nodes, 1 per blob, and vector of vector of edge indices
    net = startree(taxa) # root numbered -2
    ni = Ref(-3) # internal nodes: numbered -3,-4 etc., as done by readnewick
    ei = Ref(N+1)
    for bpart in blobparts
        weight = freq(bpart)/nnets
        # fixit: keep blob node and blob edges
        add_blobnode!(net, ni, ei, bpart.partition, weight)
    end
    for bpart in biparts
        bpart.split[N] && error("bipartition side that contains the last taxon")
        weight = freq(bpart)/nnets
        add_clusteredge_weight!(net, ni, ei, bpart.split, weight, supportaslength)
    end
    # fixit: return blob map
    return net, nothing
end

"""
    add_blobnode!(tree, ni, ei, blobpartition, weight)

Add nodes (numbered `ni` etc.) and edges (numbered `ei` etc.) in `tree`
to add `blobpartition`, assumed compatible with edges already in `tree`.
The node index `ni` is decremented, and the edge index `ei` is incremented.
Output: newly created node whose removal disconnects the taxon set into
the input blob partition.

If the blob has P ≥ 3 taxon blocks, these blocks are assumed to form a
partition of the full taxon set. If it is blob-compatible with all blobs
already in `tree`, then k nodes and k edges are added, where k is the number
of non-trivial taxon blocks.
If the blob is tree-compatible but not blob-compatible with the `tree`,
then fewer nodes and edges are added with a warning, or an error is thrown.

The blob's weight is stored in the corresponding node's `.fvalue`.

As blob partitions are agnostic about the root, the output tree should be
considered unrooted. The added edges correspond to using the last taxon as
outgroup: an edge's cluster of descendants does not contain the last taxon.

Assumptions (none are checked): `tree` is a tree, P ≥ 3, taxon blocks are
non-empty and do form a partition, and
assumptions in [`add_clusteredge!`](@ref).
"""
function add_blobnode!(
    net::PN.HybridNetwork,
    ni::Base.RefValue{Int},
    ei::Base.RefValue{Int},
    blobpartition::NTuple{P,NTuple{N,Bool}},
    weight::Number,
) where {P,N} # 2026-01: Aqua has a problem with this. "Unbound type parameters"
    # 1. create (or find) blob node: its clade is the complement of the
    #    taxon block containing the outgroup (last taxon N)
    outcluster_j = findfirst(v -> v[N], blobpartition)
    if sum(blobpartition[outcluster_j]) == 1 # trivial cluster
        lca = getroot(net)
        blobnode = lca
        # e_outj = net.edge[N] because net initialized with startree()
    else # non-trivial because P ≥ 3
        outcluster = .!blobpartition[outcluster_j]
        lca, ci = _lca_newcluster(net, outcluster)
        length(ci) == 1 && @warn("outgroup clade already in tree (or incompatible?): $outcluster")
        if length(ci) > 1
            e_outj = _resolveclade_belowlca(net, lca, ci, ni, ei, -1, false)
            blobnode = getchild(e_outj)
        end
    end
    blobnode.fvalue = weight
    # 2. create a child node below blobnode for each taxon block
    for j in 1:P
        j == outcluster_j && continue
        v = blobpartition[j]
        sum(v) == 1 && continue # skip trivial blocks of a single taxon
        node2clade_intersection_initialize(net, v)  # update .booln{2,3} for v
        node2clade_intersection_update(getroot(net))
        blobnode.booln2  || error("hmm, blob node has no taxa from the clade")
        !blobnode.booln3 || error("hmm, blob node has no taxa outside the clade")
        ci = _lca_newcluster_children(blobnode)
        length(ci) == 1 && @warn("taxon block already in tree (or incompatible?): $v")
        if length(ci) > 1
            e_j = _resolveclade_belowlca(net, blobnode, ci, ni, ei, -1, false)
        end
    end
    # fixit: also return the vector of the e_j indices from ei's.
    return blobnode
end

"""
    add_clusteredge_weight!

Find or create an edge `e` whose descendant taxa is the input cluster,
and stores the input weight in the edge's field `.y`, and as its length
with option `supportaslength=true`. Output: edge `e`.

The input cluster may already exists in the input tree.
Unlike [`add_clusteredge!`](@ref), no warning is thrown. The appropriate
edge is found, its fields (`.y` and perhaps length) are modified to store the
weight, and this edge is returned.
"""
function add_clusteredge_weight!(
    net::PN.HybridNetwork,
    ni::Base.RefValue{Int},
    ei::Base.RefValue{Int},
    bv::SplitTuple,
    weight::Number,
    supportaslength::Bool,
)
    lca, ci = _lca_newcluster(net, bv)
    length(ci) > 0 ||
        error("would not connect the new clade node to any children")
    if length(ci) == 1
        e = lca.edge[ci[1]]
        e.y = weight
        if supportaslength
            e.length = weight
        end
    else
        e = _resolveclade_belowlca(net, lca, ci, ni, ei, weight, supportaslength)
    end
    return e
end

function count_hybridclusters(blobvec::Vector{BlobFreq{N}}) where N
    hybdict = Dict{NTuple{N,Bool},Int}()
    for bf in blobvec
        for (hi, hf) in bf.hybrid
            hcluster = bf.partition[hi]
            hybdict[hcluster] = get(hybdict, hcluster, 0) + hf
        end
    end
    return hybdict
end

function _hybridsupport_cladesintree(
    hybdict::Dict{NTuple{N,Bool},Int},
    tre::PN.HybridNetwork,
    taxa::AbstractVector{<:String},
    nnets::Number,
) where {N}
    hwm = hardwiredclusters(tre, taxa) # but excludes external edges
    # map: (edge_number, reverse_direction?) => hybrid proportion
    edge2hyb = Dict{Tuple{Int,Bool},Float64}()
    for row in axes(hwm,1)
        clade = ntuple(j -> Bool(hwm[row, j+1]), N)
        if haskey(hybdict, clade)
            edge2hyb[(hwm[row,1],true)]  = hybdict[clade]/nnets
        end
        cladecomp = .!clade
        if haskey(hybdict, cladecomp)
            edge2hyb[(hwm[row,1],false)] = hybdict[cladecomp]/nnets
        end
    end
    for j in 1:N # external edges: assume numbered 1:N, from startree()
        clade = ntuple(isequal(j), N)
        if haskey(hybdict, clade)
            edge2hyb[(j,true)]  = hybdict[clade]/nnets
        end
        cladecomp = .!clade
        if haskey(hybdict, cladecomp)
            edge2hyb[(j,false)] = hybdict[cladecomp]/nnets
        end
    end
    return edge2hyb
end
