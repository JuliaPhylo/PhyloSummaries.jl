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

See also: [`consensus_treeofblobs`](@ref).

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
    update_blobcircorderfrequency!(refblobs, blobvec, taxa)
    update_hybridclusterfrequency!(refblobs, hybdict)
    # match each reference bipartition against sample frequencies
    update_bipartitionfrequency!(refbps, bpvec)
    # todo: build tables. use hwmatrix, edgemap (from above)
    # to calculate the blob-edge-indices bbei's, then
    # call blobdata_*, hybriddata_onToB and bipartdata_onToB?
    return (blob_table=blob_table, bipartition_table=bp_table, taxa=taxa)
end

"""
    update_blobcircorderfrequency!(refblobs, sampleblobs, taxa)

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
    taxa::AbstractVector{<:String},
) where N
    # todo: update the circular order
    for rb in refblobs
        matchidx, _ = findmatchingblob(sampleblobs, rb.partition)
        isnothing(matchidx) && continue # all frequencies were initialized to 0
        bf = sampleblobs[matchidx]
        freq!(rf, freq(bf))
        # copy each circular order from bf.circorder into rb.circorder
        # rb.partition[k] = bf.partition[idxmap[k]] so order
        for (co_block, co_freq) in bf.circorder
            # todo: find block order in rb, corresponding to co_block in bf
            idxmap_inrb = idxmap # wrong: fixit
            add_canonical_circularorder!(rb.circorder, idxmap_inrb, co_freq)
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

