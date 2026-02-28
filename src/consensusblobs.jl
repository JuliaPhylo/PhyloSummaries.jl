
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
- For each part in the partition, this part is a \"hybrid\" for this blob if it is
  adjacent to the blob at a hybrid node.

`P` is the number of parts (taxon blocks) in the partition; 3 or more.
`N` is the number of leaves (taxa) in the network. 3 or more.
"""
struct BlobFreq{N,P}
    """blob partition: P parts, each described as a tuple of size N
    (an immutable type, so partitions can be keys in sets or dictionaries)"""
    partition::Vector{NTuple{N,Bool}}
    "frequency of the blob partition. mutable: use freq and freq! to get/set this value."
    freq::Base.RefValue{Float64}
    "frequencies of the different circular orders for blobs with this partition"
    circorder::Dict{NTuple{P,Int},Float64}
    "frequencies of the different hybrid parts for blobs with this partition"
    hybrid::Dict{Int,Float64}
end
partitionstring(vv) = join(join.(findall.(vv),","), "|")
partitionstring(vv, itr) = join(join.([findall(vv[i]) for i in itr],","), "|")
partitionstring(obj::BlobFreq) = partitionstring(obj.partition)
partitionstring_names(vv, taxa::Vector) = join(
    (join((taxa[j] for (j,b) in enumerate(v) if b), ",") for v in vv), "|")
partitionstring_names(vv, taxa::Vector, itr) = join(
    (join((taxa[j] for (j,b) in enumerate(vv[i]) if b), ",") for i in itr), "|")

freq(obj::BlobFreq) = obj.freq[]
freq!(obj::BlobFreq, n) = obj.freq[] = Float64(n)
incrementfreq!(obj::BlobFreq, x=1) = obj.freq[] += x

Base.show(io::IO, obj::BlobFreq{N,P}) where {N,P} = print(io,
    "BlobFreq on $N taxa, $P blocks " * partitionstring(obj) *
    ", frequency $(Int(round(freq(obj))))" *
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
# tested but not used so far
function partition_to_names(io::IO, partition_numeric::AbstractString, taxa::Vector)
    join(io, (join(
        (taxa[parse(Int, jstr)] for jstr in eachsplit(str, ',')), ',')
        for str in eachsplit(partition_numeric, '|')), '|')
end
partition_to_names(p::AbstractString, taxa::Vector) = sprint(partition_to_names, p, taxa)

blobnodename(i) = "_blob$i"

iscompatible(b1::SplitFreq{N}, b2::BlobFreq{N}) where N =
    splitblobcompatible(b1.split, b2.partition) # in utils.jl
iscompatible(b1::BlobFreq{N}, b2::BlobFreq{N}) where N =
    blobcompatible(b1.partition, b2.partition)
isredundantsplit(b1::SplitFreq{N}, b2::BlobFreq{N}) where N =
    isredundantsplit(b1.split, b2.partition)

"""
    consensus_treeofblobs(networks; proportion=0,
        minimumblobdegree=4, network_weights=nothing)

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
Setting `minimumblobdegree` to 3 will cause non-trivial blobs
to be considered "interesting" even if their corresponding node in N's
tree of blobs is of degree 3.

Note that a node of degree 4 or more in the network's tree of blob may
correspond to a polytomy in N: a single node incident to m cut-edges, but
without any reticulation. These blobs are considered "non-interesting".
A cut-edge incident to such a polytomy is then non-redundant,
if the other blob it connects to is also non-interesting.

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
    suppressinfo::Bool=false,
    netweight::Union{Nothing,AbstractVector}=nothing,
)
    isempty(networks) &&
        throw(ArgumentError("No input networks: cannot get a consensus"))
    minimumblobdegree ≥ 3 ||
        throw(ArgumentError("minimumblobdegree should be 3 or higher, not $(minimumblobdegree)"))
    taxa = sort(tiplabels(networks[1]))
    nnets = length(networks)
    if !isnothing(netweight)
      length(netweight) == nnets ||
          error("there should be $nnets network weights, got $(length(netweight))")
      all(netweight .>= 0) || error("network weights should be > 0")
      nnets = sum(netweight)
    end
    blobvec, bpvec = count_blobpartitions(networks, taxa, minimumblobdegree,
        false, netweight)
    hybdict = count_hybridclusters(blobvec) # before filtering
    filter_sort_compatible_partitions!(blobvec, bpvec, nnets, proportion)
    update_hybridclusterfrequency!(blobvec, hybdict)
    tob, bbn, bbei, bpei = tree_from_blobpartitions(taxa, blobvec, bpvec, nnets, supportaslength)
    bdat, odat = blobdata_onToB(blobvec, bbn, nnets, taxa)
    hdat = hybriddata_onToB(blobvec, bbn, bbei, nnets, tob.edge, taxa)
    sdat = bipartdata_onToB(bpvec, bpei, nnets, tob.edge, taxa)
    suppressinfo || @info """
    Node & edge numbers in tables correspond to current numbers in net.
    If the network is modified or re-read from file, restore matching
    numbers with `resetnodenumbers_fromnames!` and `edgenumber`.
    """
    return (tob=tob, blob_table=bdat, circorder_table=odat,
        hybrid_table=hdat, bipartition_table=sdat, taxa=taxa)
end

"""
    consensus_level1network(networks; proportion=0,
        minimumblobdegree=4, outgroup=nothing, network_weights=nothing)

Consensus network summarizing a list of level-1 networks, by these steps:
1. A consensus tree of blobs is built as in [`consensus_treeofblobs`](@ref),
   with one node for each blob present in a majority (or in more than the
   required `proportion`) of input networks, if compatible with blobs
   of higher support, with non-redundant bipartitions supported by more
   than the `proportion` of input networks.
2. Each blob is resolved as a cycle, one after another, from highest to lowest
   supported blobs.
3. To resolve a blob, its taxon blocks are placed around a cycle in the
   circular order most-frequently found in input networks.
3. To orient the edges in the cycle, the node chosen to be hybrid is the one
   whose descendant clade has the highest support as being a hybrid clade, among
   the placements that are admissible (do not conflict with the direction of
   edges from hybrids in other cycles).
   If an `outgroup` is provided, the hybrid node is chosen among those that
   do not conflict with this outgroup taxon being an outgroup: below the root,
   and not affected (below) any reticulation.

See [`consensus_level1network_save`](@ref) to save the output.
"""
function consensus_level1network(
    networks::AbstractVector{PN.HybridNetwork};
    proportion::Number=0,
    minimumblobdegree::Int=4,
    outgroup::Union{Nothing,String}=nothing,
    suppressinfo::Bool=false,
    netweight::Union{Nothing,AbstractVector}=nothing,
)
    isempty(networks) &&
        throw(ArgumentError("No input networks: cannot get a consensus"))
    minimumblobdegree ≥ 3 ||
        throw(ArgumentError("minimumblobdegree should be 3 or higher, not $(minimumblobdegree)"))
    taxa = sort(tiplabels(networks[1]))
    isnothing(outgroup) || outgroup ∈ taxa || # early problem detection
        error("outgroup $outgroup is not in the taxon list: $taxa")
    nnets = length(networks)
    if !isnothing(netweight)
      length(netweight) == nnets ||
          error("there should be $nnets network weights, got $(length(netweight))")
      all(netweight .>= 0) || error("network weights should be > 0")
      nnets = sum(netweight)
    end
    blobvec, bpvec = count_blobpartitions(networks, taxa, minimumblobdegree,
        true, netweight)
    hybdict = count_hybridclusters(blobvec) # before filtering blobs out
    filter_sort_compatible_partitions!(blobvec, bpvec, nnets, proportion)
    update_hybridclusterfrequency!(blobvec, hybdict)
    net, bbn, bbei, bpei = tree_from_blobpartitions(taxa, blobvec, bpvec, nnets, false)
    res = expand_blobcycles!(net, bbn, bbei, blobvec, nnets, outgroup)
    bdat = blobdata_onL1(blobvec, bbn, res..., nnets, taxa)
    hdat = hybriddata_onToB(blobvec, bbn, bbei, nnets, net.edge, taxa)
    sdat = bipartdata_onToB(bpvec, bpei, nnets, net.edge, taxa)
    suppressinfo || @info """
    Node & edge numbers in tables correspond to current numbers in net.
    If the network is modified or re-read from file, restore matching
    numbers with `resetnodenumbers_fromnames!` and `edgenumber`.
    """
    return (net=net, blob_table=bdat,
        hybrid_table=hdat, bipartition_table=sdat, taxa=taxa)
end

"""
    consensus_level1network_save(result_object, rootname=nothing)

Write to files the results of [`consensus_level1network`](@ref) to 4 files,
with file names starting with `rootname` and ending with "_net.nwk" for the
consensus network (with bipartition support as non-redundant with blobs),
"_blob.csv" for the table about blobs and their support,
"_hybrid.csv" for the table about hybrid clusters and their support, and
"_bipartition.csv" for the table about bipartitions and their support
(as non-redundant with blobs).

For the hybrid and bipartition tables, the column for the edge number is not saved
because re-reading the network from the newick file would most likely lead to
different internal edge numbers, and using obsolete edge numbers could then
cause unintended errors. Node numbers will also be different when re-reading
the network from the newick file, but the original node numbers can be recovered
with [`resetnodenumbers_fromnames!`](@ref)

!!! warn "Files will be overwritten, if they already exist"

"""
function consensus_level1network_save(
    res::NamedTuple,
    root::Union{AbstractString, Nothing}=nothing
)
    isnothing(root) && error("Please provide a root for the file names")
    writenewick(res[:net], root * "_net.nwk", support=:y)
    CSV.write(root * "_blob.csv", res[:blob_table])
    CSV.write(root * "_hybrid.csv", res[:hybrid_table][ # do not save :edge column
        [:blob, :node_from, :node_to, :support_hybrid, :cluster_num, :cluster]
    ])
    CSV.write(root * "_bipartition.csv", res[:bipartition_table][ # no edge number
        [:node1, :node2, :support_nonredundant, :cluster_num, :cluster]
    ])
    return nothing
end

"""
    count_blobpartitions(networks, taxa, minimumblobdegree,
        require_level1=false, netweight=nothing)

`(blob_vec, bipart_vec)` where `blob_vec` is a vector of
[`BlobFreq{ntax}`](@ref) object (`ntax` being the number of taxa),
`bipart_vec` is a vector of [`SplitFreq{ntax}`](@ref) objects.
All input networks must have the same set of `taxa`.
If `require_level1` is true, an error is thrown the first network with a blob
whose level is not ≤ 1.

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
    require_level1::Bool=false,
    netweight::Union{Nothing,AbstractVector}=nothing,
)
    ntaxa = length(taxa)
    all(n.numtaxa == ntaxa for n in networks) ||
        throw(ArgumentError("input networks have different numbers of taxa"))
    # hardwiredclusters will error if different taxon sets
    blobvec = BlobFreq{ntaxa}[]
    bpvec = SplitFreq{ntaxa}[]  # bipartitions, frequency: if non-redundant
    nweight = (isnothing(netweight) ? i -> 1.0 : i -> Float64(netweight[i]))
    for (i,net) in enumerate(networks)
        count_blobpartitions!(blobvec, bpvec, net, taxa, minBdegree,
            require_level1, nweight(i))
    end
    return blobvec, bpvec
end

"""
    count_blobpartitions!(blobs, biparts, net, netweight, taxa, minBdegree,
        require_level1, netweight)

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
    require_level1::Bool,
    netweight::Float64,
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
            net, taxaindex, minBdegree, bc, bidx, hwmatrix, edgemap,
            require_level1, netweight)
    end
    # gather non-redundant cut-edges: from trivial bicomponents
    count_nonredundantbipartitions!(bpvec, blobdegree,
            net, taxaindex, minBdegree, hwmatrix, edgemap, netweight)
    return nothing
end

"""
    count_blobpartitions!(blobs, visitedbcc, net, taxaindex, minBdegree,
        blob, bidx, hwmatrix, edgemap, require_level1, netweight)

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
    require_level1::Bool,
    netweight::Float64,
) where N
    blobdegree = (bidx==1 ? Ref(0) : Ref(1)) # parent edge: 0 if root, 1 ow
    splits, hybrids, islevel1 = blobtaxonsetpartition!(visitedbcc, blobdegree,
        net, blob, bidx, edgemap, hwmatrix, taxaindex, N)
    if blobdegree[] < minBdegree
        return blobdegree[]
    end
    !islevel1 && require_level1 &&
        error("network with a blob of level > 1")
    isempty(splits)  && error("non-trivial blob without any split.")
    isempty(hybrids) && error("non-trivial blob without any hybrid.")
    nparts = length(splits)
    # check if this partition already exists in blobvec
    matchidx, idxmap = findmatchingblob(blobvec, splits)
    if isnothing(matchidx) # add new blob to blobvec
        if islevel1
            defaultorder = ntuple(identity, nparts)
            cofreq = Dict(defaultorder => netweight)
        else
            # level > 1: do not store circular order
            cofreq = Dict{NTuple{nparts,Int},Float64}()
        end
        hybmap = Dict(hybrids[1] => 1)
        newblob = BlobFreq{N,nparts}(splits, Ref(netweight), cofreq, hybmap)
        push!(blobvec, newblob)
    else # existing blob: increment frequencies of canonical partition slots
        bf = blobvec[matchidx]
        for k in hybrids
            canonk = idxmap[k]
            bf.hybrid[canonk] = get(bf.hybrid, canonk, 0.0) + netweight
        end
        # calculate & store circular order if level-1 blob with 1 lowest hybrid
        # a level-1 blob of >1 bicomponents could have >1 lowest hybrids
        if islevel1 && length(hybrids) == 1
            add_canonical_circularorder!(bf.circorder, idxmap, netweight)
        end
        incrementfreq!(bf, netweight)
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
    add_canonical_circularorder!(circularorder_dictionary, indexmap, netweight)

1. Find the clockwise (and counter-clockwise if necessary) circular permutation(s)
   of `indexmap`, starting at the index for value 1. For example, for vector
   `[5, 1, 3, 2, 4]`, these are: `1,3,2,4,5` and `1,5,4,2,3`.
   These are the 2 canonical ways of coding their shared circular order:
   starting from value 1 and circling in one or the other direction.
2. Add this circular order to the input dictionary: if already present
   (in clockwise or counterwise direction) then its frequency is increased by
   `netweight`.
   Otherwise, a new entry with frequency `netweight` is added to the dictionary.
"""
function add_canonical_circularorder!(
    codict::Dict{NTuple{P,Int},T},
    idxmap::AbstractVector{Int},
    netweight::T,
) where {P,T}
    startidx = findfirst(==(1), idxmap)
    isnothing(startidx) &&
        error("block 1 not found in idxmap: $idxmap (should span 1:$P)")
    itr = Iterators.flatten((startidx:P, 1:(startidx-1)))
    circorderkey  = Tuple(Iterators.map(i -> idxmap[i], itr))
    length(circorderkey) == P ||
        error("circular order key of length $(length(circorderkey)) instead of $P")
    if haskey(codict, circorderkey)
        codict[circorderkey] += netweight
    else
        itr = Iterators.flatten((startidx:-1:1, P:-1:(startidx+1)))
        reversekey = Tuple(Iterators.map(i -> idxmap[i], itr))
        if haskey(codict, reversekey)
            codict[reversekey] += netweight
        else
            codict[circorderkey] = netweight
        end
    end
end

"""
    blobtaxonsetpartition!(visitedbcc, blobdegree, net, bicomponent, bidx,
        edgemap, hardwiredclustermatrix, taxon2index, ntaxa)

Depth-first traversal of a blob `B` starting from its top bicomponent,
to collect its taxon blocks, which form a partition of the full set of `N` taxa.
Each taxon block (or split) corresponds
to a cut-edge `uv` adjacent to the blob, with `u ∈ B` in the blob and `v ∉ B`.
The taxon block for `uv` is reprented by an `N`-tuple of 0/1 values with 1 at
index `i = taxaindex[label]` if the taxon named `label` is a descendant of
`uv`, and 0 otherwise.
If `u ∈ B` is a hybrid node incident to some exit edge `uv`, then the taxon
block associated with this edge is considered "hybrid" for this blob.

Output: `(splits, hybrids, islevel1)` where
- `splits` is the blob's partition as a tuple of N-tuples, listed in a "half"
  circular order if the blob is level-1: from highest to lowest along one side
  then along the other side.
- `hybrids` is the vector of indices in `splits`, of taxon blocks
  that are hybrid for the blob.
- `islevel1` is true all bicomponent in the blob have a reticulation number ≤ 1;
  false otherwise. For binary networks, a non-trivial blob has a single
  biconnected component, and its level is its number of hybrid nodes.

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
    islevel1 = Ref(PN.getlevel(bicomp) <= 1)
    # traverse the blob starting at entrynode of the biconnected component
    blobtaxonsetpartition!(splits, hybrids, visitedbcc, blobdegree, islevel1,
        entrynode, bidx, edgemap, hwmatrix, taxaindex, net)
    nsplits = length(splits)
    if islevel1[] && length(hybrids) == 1 && hybrids[1] < nsplits
        reverse!(view(splits, (hybrids[1]+1):nsplits)) # to get circular order
    end
    if bidx > 1 # entry ≠ root: add split from the unique entry cut-edge
        outgroupsplit = splitcomplement(splits)
        if any(outgroupsplit) # empty if entry is at or above LSA (≠ root)
            # first to get circular order later: treat a node before its children
            pushfirst!(splits, outgroupsplit)
            hybrids .+= 1
        end
    end
    return splits, hybrids, islevel1[]
end

"""
    blobtaxonsetpartition!(splits, hybrids, visitedbcc, blobdegree, islevel1,
        entrynode, bidx, edgemap, hardwiredclustermatrix, taxon2index, net)

Helper to accumulate entries in `splits` and `hybrids`.
Internal fields:
- `.inte1` used but not modified: should store the index of an edge's
  biconnected component.
- `.intn1` updated to store the blob index that a visited node is in:
  index of the first biconnected component in this blob.
- `.boole1` set to `false` for all partners of an edge that will be traverered.

In `splits`, only the *children* taxon blocks are gathered, from exit cut-edges,
but recursively across all bicomponents in the blob (which may occur in
non-binary networks).
"""
function blobtaxonsetpartition!(
    splits::Vector{NTuple{N,Bool}},
    hybrids::Vector{Int},
    visitedbcc::Set{Int},
    blobdegree::Base.RefValue{Int},
    islevel1::Base.RefValue{Bool},
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
            blobtaxonsetpartition!(splits, hybrids, visitedbcc, blobdegree, islevel1,
                cn, bidx, edgemap, hwmatrix, taxaindex, net)
        elseif atotherBC && !PN.istrivial(bcc[ei]) && !(ei ∈ visitedbcc)
            push!(visitedbcc, ei)
            otherBCentry = net.vec_node[PN.entrynode_preindex(bcc[ei])]
            otherBCentry.intn1 = node.intn1
            # check if other BCC is level > 1: if so, blob is not level-1
            if PN.getlevel(bcc[ei]) > 1
                islevel1[] = false
            end
            blobtaxonsetpartition!(splits, hybrids, visitedbcc, blobdegree, islevel1,
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
    bpvec::Vector{SplitFreq{N}},
    blobdegree::Vector{Int},
    net::PN.HybridNetwork,
    taxaindex::Dict{String,Int},
    minBdegree::Int,
    hwmatrix::AbstractMatrix,
    edgemap::Dict{<:Integer,<:Integer},
    netweight::Float64,
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
        add_split!(bpvec, split, netweight)
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
            add_split!(bpvec, split, netweight)
        end
    end
    return nothing
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
    if proportion ≈ 1 # strict consensus
        filter!(v -> freq(v) ≈ nnets, blobparts)
        filter!(v -> freq(v) ≈ nnets, biparts)
    elseif threshold2 > 0 # 0 for greedy consensus: frequent case
        filter!(v -> freq(v) > threshold2, blobparts)
        filter!(v -> freq(v) > threshold2, biparts)
    end
    sort!(blobparts, by=freq, rev=false)
    sort!(biparts,   by=freq, rev=false)
    if proportion >= 0.5 # then necessarily blob-compatible, no need to check
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

Output: `(tree, blobnode_vec, blockedgeindices_vec, bipartedgeindices_vec)`
where the last 3 components list the nodes (in the tree) corresponding to each
input blob, the edges corresponding to each taxon block in each input blob, and
the edges corresponding to each non-redundant bipartition.

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
    blobnode = PN.Node[]
    blobedges = Vector{Int}[] # edge indices, one for each taxon block
    net = startree(taxa) # leaves numbered 1:N and root numbered N+1
    ni = Ref(N+2)
    ei = Ref(N+1)
    for (bnum, bpart) in enumerate(blobparts)
        weight = freq(bpart)/nnets
        bn, be = add_blobnode!(net, ni, ei, bpart.partition, weight)
        bn.name *= blobnodename(bnum) # append: keep info about node number
        push!(blobnode, bn)
        push!(blobedges, be)
    end
    bpedges = Int[] # edge indices, one for each bipartition
    for bpart in biparts
        bpart.split[N] && error("bipartition side that contains the last taxon")
        weight = freq(bpart)/nnets
        e = add_clusteredge_weight!(net, ni, ei, bpart.split, weight, supportaslength)
        push!(bpedges, e.number)
    end
    return net, blobnode, blobedges, bpedges
end

"""
    add_blobnode!(tree, ni, ei, blobpartition, weight)

Add nodes (numbered `ni` etc.) and edges (numbered `ei` etc.) in `tree`
to add `blobpartition`, assumed compatible with edges already in `tree`.
The node index `ni` is decremented, and the edge index `ei` is incremented.

Output: `(n, ei)` where
- `n` is the newly created node in `tree` whose removal disconnects the
  taxon set into the input blob partition, and
- `ei` is a vector of edge indices in `tree.edge`: `tree.edge[e[j]]` is
  the edge whose removal disconnects taxon block `j` of `blobpartition`
  from its complement.

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
non-empty and do form a partition, assumptions in [`add_clusteredge!`](@ref),
and taxon `j` is the node number `j` in `tree` incident to `tree.edge[j]`,
as built from [`startree`](@ref).
"""
function add_blobnode!(
    net::PN.HybridNetwork,
    ni::Base.RefValue{Int},
    ei::Base.RefValue{Int},
    blobpartition::Vector{NTuple{N,Bool}},
    weight::Number,
) where N
    P = length(blobpartition)
    bei = Vector{Int}(undef, P) # blob edge indices
    # 1. create (or find) blob node: its clade is the complement of the
    #    taxon block containing the outgroup (last taxon N)
    outcluster_j = findfirst(v -> v[N], blobpartition)
    if sum(blobpartition[outcluster_j]) == 1 # trivial cluster
        lca = getroot(net)
        blobnode = lca
        bei[outcluster_j] = N # because net initialized with startree()
        # e_outj = net.edge[N]
    else # non-trivial because P ≥ 3
        outcluster = .!blobpartition[outcluster_j]
        lca, ci = _lca_newcluster(net, outcluster)
        if length(ci) == 1 # edge may already exist if implied (redundant) by 2 blobs
            e_outj = lca.edge[ci[1]]
        else # length(ci) > 1
            e_outj = _resolveclade_belowlca(net, lca, ci, ni, ei, -1, false)
        end
        blobnode = getchild(e_outj)
        bei[outcluster_j] = e_outj.number # because net.edge[k].number = k
    end
    blobnode.fvalue = weight
    # 2. create a child node below blobnode for each taxon block
    for j in 1:P
        j == outcluster_j && continue
        v = blobpartition[j]
        if sum(v) == 1 # skip trivial blocks of a single taxon
            bei[j] = findfirst(v) # taxon i: has edge index & number i, from startree
            continue
        end
        node2clade_intersection_initialize(net, v)  # update .booln{2,3} for v
        node2clade_intersection_update(getroot(net))
        blobnode.booln2  || error("hmm, blob node has no taxa from the clade")
        !blobnode.booln3 || error("hmm, blob node has no taxa outside the clade")
        ci = _lca_newcluster_children(blobnode)
        if length(ci) == 1
            bei[j] = blobnode.edge[ci[1]].number
            # @warn "taxon block already in tree (or incompatible?): $v"
        else # length(ci) > 1
            bei[j] = ei[]
            e_j = _resolveclade_belowlca(net, blobnode, ci, ni, ei, -1, false)
        end
    end
    return blobnode, bei
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

# for each hybrid clade, sum its frequency over all blobs that may have it
function count_hybridclusters(blobvec::Vector{BlobFreq{N}}) where N
    hybdict = Dict{NTuple{N,Bool},Float64}()
    for bf in blobvec
        for (hi, hf) in bf.hybrid
            hcluster = bf.partition[hi]
            hybdict[hcluster] = get(hybdict, hcluster, 0) + hf
        end
    end
    return hybdict
end

"""
    update_hybridclusterfrequency!(blobfreq_vector, clusterfreq_dict)

In each `BlobFreq` item `blob` in `blobfreq_vector`, replace hybrid frequencies
in `blob.hybrid` by those in `clusterfreq_dict`.
Note that the list of hybrid clusters in `blob` may be expanded with clusters
not yet listed but present in `clusterfreq_dict`, for any cluster compatible
with the blob partition. Also, clusters are rooted here (hybrid descendants).
Assumption: any part of a `blob` is a cluster (key) in `clusterfreq_dict`.
"""
function update_hybridclusterfrequency!(
    blobvec::Vector{BlobFreq{N}},
    hybdict::Dict{NTuple{N,Bool},Float64},
) where N
    for bf in blobvec
        for (hi,cluster) in enumerate(bf.partition)
            if haskey(hybdict, cluster)
                bf.hybrid[hi] = hybdict[cluster] # update of add item
            end
        end
    end
    return nothing
end

function blobdata_onToB(
    blobparts::Vector{BlobFreq{N}},
    blobnode::Vector{PN.Node},
    nnets::Number,
    taxa::Vector{<:AbstractString},
) where N
    nB = length(blobparts)
    @assert nB == length(blobnode)
    bitr = ((i,blobparts[i]) for i in nB:-1:1) # from most to least frequent blob
    blob_data = (blob = [i for (i,b) in bitr],
        degree = [length(b.partition) for (i,b) in bitr],
        node = [blobnode[i].number for (i,b) in bitr],
        support_partition = [freq(b)/nnets for (i,b) in bitr],
        partition_num = [partitionstring(b) for (i,b) in bitr],
        partition = [partitionstring_names(b.partition, taxa) for (i,b) in bitr])
    itr = ((i,b,co,f) for (i,b) in bitr for (co,f) in b.circorder)
    co_data = (blob = [x[1] for x in itr],
        order = [x[3] for x in itr],
        support_circorder = [x[4]/nnets for x in itr],
        partition_num = [partitionstring(b.partition, co) for (i,b,co,f) in itr],
        partition = [partitionstring_names(b.partition, taxa, co) for (i,b,co,f) in itr])
    return blob_data, co_data
end

function blobdata_onL1( # for consensus level-1 network
    blobparts::Vector{BlobFreq{N}},
    blobnode::Vector{PN.Node},
    o_bs::Vector,
    h_bs::Vector,
    h_num::Vector,
    h_blk::Vector,
    nnets::Number,
    taxa::Vector{<:AbstractString},
) where N
    nB = length(blobparts)
    @assert nB == length(blobnode)
    bitr = ((i,blobparts[i]) for i in nB:-1:1) # from most to least frequent blob
    blob_data = (
        blob = [i for (i,b) in bitr],
        degree = [length(b.partition) for (i,b) in bitr],
        node = [blobnode[i].number for (i,b) in bitr],
        hybrid = [h_num[i] for (i,b) in bitr],
        support_partition = [freq(b)/nnets for (i,b) in bitr],
        support_circorder = [o_bs[i] for (i,b) in bitr],
        support_hybrid = [h_bs[i] for (i,b) in bitr],
        partition_num = [partitionstring(b) for (i,b) in bitr],
        hybrid_cluster_num = [splitstring(b.partition[h_blk[i]]) for (i,b) in bitr],
        partition = [partitionstring_names(b.partition, taxa) for (i,b) in bitr],
        hybrid_cluster = [splitstring_names(b.partition[h_blk[i]], taxa) for (i,b) in bitr],
    )
    return blob_data
end

function hybriddata_onToB(
    blobparts::Vector{BlobFreq{N}},
    blobnode::Vector{PN.Node},
    blobedges::Vector{Vector{Int}},
    nnets::Number,
    netedge::Vector{PN.Edge},
    taxa::Vector{<:AbstractString},
) where N
    nB = length(blobparts)
    @assert nB == length(blobedges)
    bitr = ((i,blobparts[i]) for i in nB:-1:1) # from most to least frequent blob
    itr = ((i,b,h,f) for (i,b) in bitr for (h,f) in b.hybrid)
    bnum = [x[1] for x in itr]
    hedgenum = [blobedges[i][h] for (i,b,h,f) in itr]
    nE = length(hedgenum)
    pnum = Vector{Int}(undef, nE) # parent & child node numbers: if
    cnum = Vector{Int}(undef, nE) # edge directed p-->c -> hybrid clade
    for (j,bi,ei) in zip(1:nE, bnum, hedgenum)
        nn = netedge[ei].node
        from1 = nn[1] === blobnode[bi] || nn[1].intn1 == bi # .intn1 not set for ToB
        pnum[j] = nn[(from1 ? 1 : 2)].number
        cnum[j] = nn[(from1 ? 2 : 1)].number
    end
    hybrid_data = (blob = bnum, node_from = pnum, node_to = cnum, edge = hedgenum,
        support_hybrid = [x[4]/nnets for x in itr],
        cluster_num = [splitstring(b.partition[h]) for (i,b,h,f) in itr],
        cluster = [splitstring_names(b.partition[h], taxa) for (i,b,h,f) in itr],
    )
    return hybrid_data
end

function bipartdata_onToB(
    biparts::Vector{SplitFreq{N}},
    biedges::Vector{Int},
    nnets::Number,
    netedge::Vector{PN.Edge},
    taxa::Vector{<:AbstractString},
) where N
    nB = length(biparts)
    @assert nB == length(biedges)
    bitr = ((i,biparts[i]) for i in nB:-1:1) # from most to least frequent
    enum = [biedges[i] for (i,b) in bitr]
    pnum = Vector{Int}(undef, length(enum)) # parent & child node numbers
    cnum = Vector{Int}(undef, length(enum))
    for i in nB:-1:1
        ee = netedge[enum[i]]
        pnum[i] = getparent(ee).number
        cnum[i] = getchild(ee).number
    end
    dat = (node1 = pnum, node2 = cnum, edge = enum,
        support_nonredundant = [freq(b)/nnets for (i,b) in bitr],
        cluster_num = [splitstring(b) for (i,b) in bitr],
        cluster = [splitstring_names(b, taxa) for (i,b) in bitr],)
    return dat
end

"""
    expand_blobcycles!(net, blobnode_vector, edgeblockindices_vector,
        blobpartition_vector, nnets, outgroup=nothing)

Modify each blob node in `net` into a cycle. If specified, `net` is rooted
at the outgroup taxon.
Blobs are expanded sequentially, from most frequent to least frequent.
Each blob node is expanded using a circular order of highest frequency.
The hybrid node in the cycle is one of highest frequency among those compatible
with previous expanded blobs, or some child of the blob node otherwise.

Weights (frequency / `nnets`) of blobs, circular order and hybrid clades are
stored in specific nodes' `.fvalue`:
- The blob weight should be stored in the fvalue of the blob node
  if the input `net` is originally from [`tree_from_blobpartitions`](@ref).
  The original blob node will become the entry (or root) of the blob cycle.
- The weight of the circular order is stored in the fvalue of all
  non-entry tree nodes in the cycle.
- The weight of a hybrid clade is stored the hybrid node's fvalue.

In each cycle, the edges' `.inte1` and the nodes' `intn1` are set to a
blob identifier (its index in the input vectors).
"""
function expand_blobcycles!(
    net::PN.HybridNetwork,
    blobnode::Vector{PN.Node},
    blobedges::Vector{Vector{Int}},
    blobparts::Vector{BlobFreq{N}},
    nnets::Number,
    outgroup::Union{Nothing, AbstractString},
) where N
    fixdirection = !isnothing(outgroup)
    if fixdirection
        oldoutedgei = findfirst(n -> n.name==outgroup, net.node) # assume edge i to leaf i
        rootatnode!(net, oldoutedgei; index=true) # creates degree-2 node
        newoutedgei = length(net.edge) # last=new edge: incident to blob if {outgroup}=block
        replace!.(blobedges, oldoutedgei => newoutedgei; count=1)
    end
    nnum = Ref(maximum(n.number for n in net.node)+1)
    enum = Ref(maximum(e.number for e in net.edge)+1)
    nB = length(blobnode)
    o_bs = Vector{Float64}(undef, nB) # order: bootstrap support
    h_bs = Vector{Float64}(undef, nB) # hybrid: bootstrap support
    h_num = Vector{Int}(undef, nB) # hybrid: node number
    h_blk = Vector{Int}(undef, nB) # hybrid: block number in the multipartition
    @assert nB == length(blobedges) == length(blobparts)
    bitr = ((i,blobparts[i]) for i in nB:-1:1) # from most to least frequent blob
    for (i,b) in bitr
        o_bs[i], h_bs[i], h_num[i], h_blk[i] = expand_blobcycleat!(net,
            nnum, enum, i, blobnode[i],blobedges[i],b, nnets, fixdirection)
    end
    return o_bs, h_bs, h_num, h_blk
end
function expand_blobcycleat!(
    net::PN.HybridNetwork,
    nnum::Base.RefValue{Int},
    enum::Base.RefValue{Int},
    bnum::Integer,
    bnode::PN.Node,
    bedges::Vector{Int},
    bpart::BlobFreq{N,P},
    nnets::Number,
    fixdirection::Bool,
) where {N,P}
    # 1. find a taxon block / edge to be the hybrid block
    hblock = argmax(bpart.hybrid) # most frequent hybrid block
    hedge = net.edge[bedges[hblock]]
    hbelowblob = isparentof(bnode, hedge)
    if !hbelowblob && !fixdirection && hedge.containroot
        rootatnode!(net, bnode) # re-root at the blob node
        hbelowblob = isparentof(bnode, hedge)
    end
    if !hbelowblob # then find another block to be hybrid
        priorh = hblock
        if length(bpart.hybrid) == 1 # then pick block 1, or 2 if prior was 1
            hblock = (priorh == 1 ? 2 : 1)
        else # pick second most frequent hybrid block
            hblock = argmax(k -> (k==priorh ? 0 : bpart.hybrid[k]), keys(bpart.hybrid))
        end
        hedge = net.edge[bedges[hblock]]
        isparentof(bnode, hedge) || error("blob node with 2 parents?")
    end
    hweight = get(bpart.hybrid, hblock, 0)/nnets
    # 2. find a block / edge to remain connected to blob node
    if isrootof(bnode, net) # no parent: pick block 1 (or 2)
        pblock = (hblock == 1 ? 2 : 1)
        pedge = net.edge[bedges[pblock]]
    else # parent edge (and its block) will stay incident to blob node
        pedge = getparentedge(bnode)
        pblock = findfirst(i -> net.edge[i] === pedge, bedges)
        isnothing(pblock) && error("top block above blob node not found")
    end
    # 3. create the cycle, in most-frequent circular order
    #    P-1 new nodes: one for each non-parent edge
    #    P   new edges: P-2 tree edges + 2 hybrid edges
    blockorder = argmax(bpart.circorder)
    circweight = bpart.circorder[blockorder] / nnets
    ii_p = findfirst(isequal(pblock), blockorder)
    ii_h = findfirst(isequal(hblock), blockorder)
    neighbor = nothing
    downward = ii_h < ii_p # used later to orient the new edges
    containroot = pedge.containroot
    for ii in Iterators.flatten((1:P, [1]))
        k = blockorder[ii]
        ee = net.edge[bedges[k]]
        ishyb_n = (ii==ii_h)
        newnode = bnode # good for parent block only
        if ii!=ii_p     # for non-parent blocks:
          if isnothing(neighbor) || ii!=1 # first P blocks: create new node
            newnode = PN.Node(nnum[], false, ishyb_n,
                (ishyb_n ? hweight : circweight), [ee]) # fvalue, edge
            newnode.name = (ishyb_n ? "H$(nnum[])" : "_$(nnum[])")
            newnode.intn1 = bnum
            nnum[] += 1
            PN.removeEdge!(bnode, ee)
            ee.node[findfirst(x -> x === bnode, ee.node)] = newnode
            PN.pushNode!(net, newnode)
          else # last iteration: back to block 1
            newnode = getparent(ee) # created at first iteration
          end
        end
        newnode.intn1 = bnum
        if !isnothing(neighbor) # last P blocks: create new edge
            ishyb_e = ishyb_n || neighbor.hybrid
            gam = (ishyb_e ? -1.0 : 1.0)
            ismajor = !ishyb_n
            newedge = PN.Edge(enum[], -1., ishyb_e, # length=-1
                -1.,-1., gam, [newnode,neighbor], downward, # y=-1, z=-1
                ismajor, bnum, containroot, true, false)
            enum[] += 1
            push!(neighbor.edge, newedge)
            push!(newnode.edge,  newedge)
            PN.pushEdge!(net, newedge)
        end
        neighbor = newnode
        if ishyb_n || ii==ii_p
            downward = !downward
        end
    end
    if containroot # to update edges' containroot. ischild1 already correct
        PN.traverseDirectEdges!(getparent(hedge), hedge, false)
    end
    return circweight, hweight, getparent(hedge).number, hblock
end
