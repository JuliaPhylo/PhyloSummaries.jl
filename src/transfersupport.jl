"""
    blobtransferdistance!(C, blob1, blob2)

Compute the transfer distance between two partitions `blob1` and `blob2` of the
same set of `N` taxa, using `C` as a pre-allocated `N×N` (or larger) `Int` matrix
buffer to avoid allocations. `C` is mutated.
"""
function blobtransferdistance!(
    C::AbstractMatrix{Int},
    blob1::AbstractVector{NTuple{N,Bool}},
    blob2::AbstractVector{NTuple{N,Bool}},
) where N
    k1 = length(blob1)
    k2 = length(blob2)
    (k1 == 0 || k2 == 0) && return 0
    kmax = max(k1, k2)
    Cv = view(C, 1:kmax, 1:kmax)
    fill!(Cv, 0)
    for i in 1:kmax
        ni = i <= k1 ? count(blob1[i]) : 0
        for j in 1:kmax
            if j <= k2
                nj = count(blob2[j])
                isect = i <= k1 ? count(t -> blob1[i][t] & blob2[j][t], 1:N) : 0
            else
                nj = 0
                isect = 0
            end
            Cv[i,j] = ni + nj - 2isect
        end
    end
    _, mincost = hungarian(Cv)
    return mincost
end

"""
    blobtransferdistance(blob1, blob2)

"""
function blobtransferdistance(
    blob1::AbstractVector{NTuple{N,Bool}},
    blob2::AbstractVector{NTuple{N,Bool}},
) where N
    kmax = max(length(blob1), length(blob2))
    C = Matrix{Int}(undef, kmax, kmax)
    return blobtransferdistance!(C, blob1, blob2)
end


"""
    transferindex!(C, dist_cache, ref_idx, refblob, net_blobidxs, blobvec)

Minimum transfer distance between reference blob `refblob` (at index `ref_idx`)
and the blobs in a sample network, represented as `net_blobidxs` (indices into
`blobvec`). Distances are lazily computed via `blobtransferdistance!` and stored
in `dist_cache[(ref_idx, j)]` for reuse across networks.
"""
function transferindex!(
    C::AbstractMatrix{Int},
    dist_cache::Dict{Tuple{Int,Int},Int},
    ref_idx::Int,
    refblob::AbstractVector{NTuple{N,Bool}},
    net_blobidxs::AbstractVector{Int},
    blobvec::Vector{BlobFreq{N}},
) where N
    isempty(net_blobidxs) && return N
    mindist = typemax(Int)
    for j in net_blobidxs
        key = (ref_idx, j)
        d = get(dist_cache, key, -1)
        if d < 0
            d = blobtransferdistance!(C, refblob, blobvec[j].partition)
            dist_cache[key] = d
        end
        d < mindist && (mindist = d)
        mindist == 0 && break
    end
    return mindist
end


"""
    blobtransfer_support(networks, referencenet; minimumblobdegree=3)
 
Compute transfer support for each blob partition in `referencenet` based on
its minimum transfer distance to blob partitions found in `networks`.

"""
function blobtransfer_support(
    networks::AbstractVector{PN.HybridNetwork},
    referencenet::PN.HybridNetwork;
    minimumblobdegree::Int=3,
)
    isempty(networks) &&
        throw(ArgumentError("No input networks: cannot compute transfer support"))
    taxa = sort(tiplabels(referencenet))
    ntaxa = length(taxa)

    refblobs = BlobFreq{ntaxa}[]
    count_blobpartitions!(refblobs, SplitFreq{ntaxa}[], referencenet,
        taxa, minimumblobdegree, false, 0.0)
    nref = length(refblobs)
    nnet = length(networks)

    # build global deduplicated blobvec with per-network blob indices
    net_blobidx = Vector{Int}[]
    blobvec, _ = count_blobpartitions(networks, taxa,
        minimumblobdegree; net_blobidx)

    # TODO: benchmark Dict vs Matrix cache approach
    C = Matrix{Int}(undef, ntaxa, ntaxa)
    dist_cache = Dict{Tuple{Int,Int},Int}()

    ti_sum = zeros(Float64, nref)
    for idxs in net_blobidx
        for r in 1:nref
            ti_sum[r] += transferindex!(C, dist_cache, r,
                refblobs[r].partition, idxs, blobvec)
        end
    end
    return (refblobs=refblobs, transfer_index_sum=ti_sum,
        nnet=nnet, taxa=taxa)
end