"""
    blobpartitions_support(networks, referencenet;
        minimumblobdegree=4, network_weights=nothing)

Calculate the support for blob partitions and related features (circular orders,
hybrid clades, bipartitions non-redundant with a blob) for the blobs present
in a reference network `referencenet`,
based on their frequency in a sample of `networks`.

Sample networks are weighted equally by default, unless a vector of
`network_weights` is provided (of same length as the vector of sample networks).

Output: a `NamedTuple` of tables, named as follows
- `:blob_table`: support for each blob partition in the reference network,
  including support for the circular order of the taxon blocks if the sample
  networks are of level 1.
- `:hybrid_table`: support for each taxon block to be below a blob's
  lowest hybrid node in a sample network.
- `:bipartition_table`: support for each non-redundant bipartition (cut-edge)
  in the reference network.
- `:taxa`: sorted taxon labels shared by all networks. `taxa[i]` is indicated
  as taxon `i` in the other tables.

See also:
[`PhyloNetworks.treeedges_support`](@extref),
[`PhyloNetworks.hybridclades_support`](@extref),
[`consensus_treeofblobs`](@ref),
[`PhyloSummaries.count_blobpartitions`](@ref).

# example

```jldoctest
todo
```
"""
function blobpartitions_support(
    networks::AbstractVector{PN.HybridNetwork},
    referencenet::PN.HybridNetwork;
    minimumblobdegree::Int=4,
    netweight::Union{Nothing,AbstractVector}=nothing,
)
    isempty(networks) &&
        throw(ArgumentError("No input networks: cannot compute reference support"))
    minimumblobdegree ≥ 3 ||
        throw(ArgumentError("minimumblobdegree should be 3 or higher, not $minimumblobdegree"))
    taxa = sort(tiplabels(referencenet))
    ntaxa = length(taxa)
    nnets = length(networks)
    if !isnothing(netweight)
      length(netweight) == nnets ||
          error("there should be $nnets network weights, got $(length(netweight))")
      all(netweight .>= 0) || error("network weights should be > 0")
      nnets = sum(netweight)
    end
    blobvec, bpvec = count_blobpartitions(networks, taxa, minimumblobdegree)
    hybdict = count_hybridclusters(blobvec)
    refblobs = BlobFreq{ntaxa}[]
    refbps = SplitFreq{ntaxa}[]
    hwmatrix, edgemap = count_blobpartitions!(refblobs, refbps, referencenet,
        taxa, minimumblobdegree, false, 0.0) # frequencies initialized at 0
    update_blobcircorderfrequency!(refblobs, blobvec)
    update_hybridclusterfrequency!(refblobs, hybdict)
    update_bipartitionfrequency!(refbps, bpvec)
    bbn, bbei_nums = _blobnode_blobedges(referencenet, hwmatrix, edgemap, refblobs, taxa)
    # convert edge numbers to edge indices for the existing table builders
    edgenum2idx = Dict(e.number => i for (i, e) in pairs(referencenet.edge))
    bbei = [Int[edgenum2idx[n] for n in nums] for nums in bbei_nums]
    bpei = _bipartition_edgeindices(refbps, referencenet, hwmatrix, edgemap, edgenum2idx, taxa)
    # build tables using existing builders
    bdat, odat = blobdata_onToB(refblobs, bbn, nnets, taxa)
    hdat = hybriddata_onToB(refblobs, bbn, bbei, nnets, referencenet.edge, taxa)
    sdat = bipartdata_onToB(refbps, bpei, nnets, referencenet.edge, taxa)
    return (blob_table=bdat, circorder_table=odat,
        hybrid_table=hdat, bipartition_table=sdat, taxa=taxa)
end

"""
    _bipartition_edgeindices(refbps, net, hwmatrix, edgemap, edgenum2idx, taxa)

For each non-redundant bipartition in `refbps`, find the corresponding edge
index in `net.edge` by matching the split against the `hwmatrix`.
"""
function _bipartition_edgeindices(
    refbps::Vector{SplitFreq{N}},
    net::PN.HybridNetwork,
    hwmatrix::AbstractMatrix,
    edgemap::Dict{<:Integer,<:Integer},
    edgenum2idx::Dict{Int,Int},
    taxa::AbstractVector{<:String},
) where N
    bpei = Int[]
    for rb in refbps
        matched = false
        csplit = .!rb.split
        for row in axes(hwmatrix, 1)
            rowsplit = cluster_fromHmatrix(hwmatrix, row, N)
            if rowsplit == rb.split || rowsplit == csplit
                push!(bpei, edgenum2idx[hwmatrix[row, 1]])
                matched = true
                break
            end
        end
        matched || error("bipartition $(rb.split) not found in hwmatrix")
    end
    return bpei
end

"""
    update_blobcircorderfrequency!(refblobs, sampleblobs)

Update the frequency and circular order dictionary of each "reference" blob
partition blob (item in `refblobs`), to match that from `sampleblobs`.
A reference blob partition not found in `sampleblobs` is not modified --
which assumes that its frequencies had been initialized to 0.

Note the blobs hybrid dictionaries are *not* updated. For this, see instead
[`update_hybridclusterfrequency!`](@ref) to use the frequencies of hybrid
clades aggregated over all sample blobs
(a hybrid clade can originate below different blobs, of different partitions).
"""
function update_blobcircorderfrequency!(
    refblobs::Vector{BlobFreq{N}},
    sampleblobs::Vector{BlobFreq{N}},
) where N
    for rb in refblobs
        matchidx, idxmap = findmatchingblob(sampleblobs, rb.partition)
        isnothing(matchidx) && continue # all frequencies were initialized to 0
        bf = sampleblobs[matchidx]
        freq!(rb, freq(bf))
        # copy circular orders from bf into rb, remapping block indices.
        # rb.partition[k] = bf.partition[idxmap[k]], so the inverse gives us:
        # bf block index j → rb block index invmap[j]
        P = nblocks(rb) # also = length(idxmap)
        invmap = Vector{Int}(undef, P)
        for k in 1:P
            invmap[idxmap[k]] = k
        end
        for (co_key, co_freq) in bf.circorder
            remapped = [invmap[j] for j in co_key]
            add_canonical_circularorder!(rb.circorder, remapped, co_freq)
        end
    end
    return nothing
end

"""
    update_bipartitionfrequency!(refbps, samplebps)

Update the frequency of each bipartition in `refbps` to match that from
`samplebps`. A bipartition not found in `samplebps` is not modified --
which assumes that its frequencies had been initialized to 0.

Bipartitions are considered unrooted: so a bipartition in `samplebps` is a
match if either its clade or its complement matches that from `refbps`.
"""
function update_bipartitionfrequency!(
    refbps::Vector{SplitFreq{N}},
    samplebps::Vector{SplitFreq{N}},
) where N
    for rb in refbps
        i = findfirst(bp -> bp.split == rb.split || bp.split == .!rb.split, samplebps)
        isnothing(i) && continue
        freq!(rb, freq(samplebps[i]))
    end
    return nothing
end

"""
    _blobnode_blobedges

Output: `(blobnode, blobedge_numbers)` where
- `blobnode` is the vector of entry nodes, one per blob and
- `blobedge_numbers` has one vector per blob, each being a vector with
  one edge number per taxon block.

**Warning**: uses `net.partition` and internal fields `.inte1` and `.intn1`.
These fields should store the edge's bicomponent number, and the node's blob
number (if non-trivial, 0 if the node is a trivial blob).
They do if `net` was already traversed by [`count_blobpartitions`](@ref),
which calls `process_biconnectedcomponents!` to build `net.partition`.
The blob number should also be its index in the input `refblobs`.
"""
function _blobnode_blobedges(
    net::PN.HybridNetwork,
    hwmatrix::AbstractMatrix,
    edgemap::Dict{<:Integer,<:Integer},
    refblobs::Vector{BlobFreq{N}},
    taxa::AbstractVector{<:String},
) where N
    blobnode = PN.Node[]
    nblobs = length(refblobs)
    blobedge_nums = [zeros(Int, nblocks(bb)) for bb in refblobs]
    intn1_to_blobidx = Dict{Int,Int}() # .intn1 (bicomp index) → position in refblobs
    i_prev = 0
    # fixit: I think it's incorrect below, generally.
    # intn1_to_blobidx can be incomplete when trying to add a cut-edge.
    # show one example where intn1 != blob ID = blob index in refblobs
    for p in net.partition
        if PN.istrivial(p) # add cut-edge to blobedge_nums?
            ee = p.edges[1]
            b1_i = ee.node[1].intn1
            b2_i = ee.node[2].intn1
            (b1_i == 0 && b2_i == 0) && continue # non-redundant bipart, could be external
            enum = ee.number
            row = get(edgemap, enum, nothing)
            if isnothing(row)
                cn = getchild(ee)
                cn.leaf || error("cut-edge without a row in hwmatrix, yet child isn't a leaf")
                i0 = findfirst(isequal(cn.name), taxa)
                split = ntuple(isequal(i0), N)
            else
                split = cluster_fromHmatrix(hwmatrix, row, N)
            end
            vsplitmatch(v) = v == split || all(v .!== split)
            for b_i in (b1_i, b2_i)
                b_i == 0 && continue # adjacent to trivial blob
                bi = get(intn1_to_blobidx, b_i, nothing)
                isnothing(bi) && continue # blob not interesting (degree < minBdegree)
                block_j = findfirst(vsplitmatch, refblobs[bi].partition)
                isnothing(block_j) && error("cut-edge $enum, split $split: no matching taxon block in $(refblobs[bi])")
                blobedge_nums[bi][block_j] = enum
            end
        else # add blobnode
            nn = net.vec_node[PN.entrynode_preindex(p)]
            # nn.intn1 ≤ i_prev if 2+ bicomps in same blob (0 if trivial blob)
            nn.intn1 <= i_prev && continue
            push!(blobnode, nn)
            intn1_to_blobidx[nn.intn1] = length(blobnode)
            i_prev = nn.intn1
        end
    end
    # fill in outgroup/complement blocks: entry cut-edge above each blob
    for (bi, bn) in pairs(blobnode)
        any(blobedge_nums[bi] .== 0) || continue
        # find the parent edge leading outside the blob
        for e in bn.edge
            if PN.ischildof(bn, e)
                j = findfirst(blobedge_nums[bi] .== 0)
                blobedge_nums[bi][j] = e.number
                break
            end
        end
    end
    any(any(enums .== 0) for enums in blobedge_nums) &&
        error("some blob taxon blocks have no matching edge")
    length(blobnode) == nblobs ||
        error("found $(length(blobnode)) blobs instead of $nblobs")
    return blobnode, blobedge_nums
end
