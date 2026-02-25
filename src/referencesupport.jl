"""
    referencesupport(networks, referencenet; minimumblobdegree=4)

Calculate the support for features (blob partitions, bipartitions, circular
orders, hybrid clades) present in a reference network `referencenet`, based on
their frequency in a sample of `networks`.


# Output

A `NamedTuple` with:
- `blob_table`: support for each blob partition in the reference network.
- `bipartition_table`: support for each non-redundant bipartition (cut-edge)
  in the reference network.
- `taxa`: sorted taxon labels shared by all networks.

See also: [`consensus_treeofblobs`](@ref), [`count_blobpartitions`](@ref)
"""
function referencesupport(
    networks::AbstractVector{PN.HybridNetwork},
    referencenet::PN.HybridNetwork;
    minimumblobdegree::Int=4,
)
    isempty(networks) &&
        throw(ArgumentError("No input networks: cannot compute reference support"))
    minimumblobdegree ≥ 3 ||
        throw(ArgumentError("minimumblobdegree should be 3 or higher, not $minimumblobdegree"))
    taxa = sort(tiplabels(referencenet))

    nnets = length(networks)
    blobvec, bpvec = count_blobpartitions(networks, taxa, minimumblobdegree)
    refblobs, refbps = count_blobpartitions([referencenet], taxa, minimumblobdegree)
    # match each reference blob partition against sample frequencies
    blob_table = _refsupport_blobs(refblobs, blobvec, nnets, taxa)
    # match each reference bipartition against sample frequencies
    bp_table = _refsupport_bipartitions(refbps, bpvec, nnets, taxa)
    return (blob_table=blob_table, bipartition_table=bp_table, taxa=taxa)
end

"""
    _refsupport_blobs(refblobs, sampleblobs, nnets, taxa)

For each blob partition in the reference network, find the matching partition
in `sampleblobs` and report its frequency and support (frequency / nnets).
"""
function _refsupport_blobs(
    refblobs::Vector{BlobFreq{N}},
    sampleblobs::Vector{BlobFreq{N}},
    nnets::Int,
    taxa::AbstractVector{<:String},
) where N
    partitions  = String[]
    nparts      = Int[]
    frequencies = Int[]
    supports    = Float64[]
    for rb in refblobs
        matchidx, _ = findmatchingblob(sampleblobs, rb.partition)
        f = isnothing(matchidx) ? 0 : freq(sampleblobs[matchidx])
        push!(partitions, partitionstring_names(rb.partition, taxa))
        push!(nparts, length(rb.partition))
        push!(frequencies, f)
        push!(supports, f / nnets)
    end
    return (partition=partitions, degree=nparts, frequency=frequencies,
            support=supports)
end

"""
    _refsupport_bipartitions(refbps, samplebps, nnets, taxa)

For each non-redundant bipartition in the reference network, find the matching
bipartition in `samplebps` (unrooted: either orientation) and report its
frequency and support (frequency / nnets).
"""
function _refsupport_bipartitions(
    refbps::Vector{SplitFreq{N}},
    samplebps::Vector{SplitFreq{N}},
    nnets::Int,
    taxa::AbstractVector{<:String},
) where N
    clusters    = String[]
    frequencies = Int[]
    supports    = Float64[]
    for rb in refbps
        i = findfirst(bp -> bp.split == rb.split || bp.split == .!rb.split, samplebps)
        f = isnothing(i) ? 0 : freq(samplebps[i])
        push!(clusters, splitstring_names(rb.split, taxa))
        push!(frequencies, f)
        push!(supports, f / nnets)
    end
    return (cluster=clusters, frequency=frequencies, support=supports)
end

