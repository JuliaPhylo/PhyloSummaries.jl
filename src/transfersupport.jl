"""
    blobtransferdistance!(C, blob_a, blob_b)

Compute the transfer distance between two partitions `blob_a` and `blob_b` of the
same set of `N` taxa, using `C` as a pre-allocated `N×N` (or larger) `Int` matrix
buffer to avoid allocations. `C` is mutated.
"""
function blobtransferdistance!(
    C::AbstractMatrix{Int},
    blob_a::AbstractVector{NTuple{N,Bool}},
    blob_b::AbstractVector{NTuple{N,Bool}},
) where N
    k1 = length(blob_a)
    k2 = length(blob_b)
    (k1 == 0 || k2 == 0) && return 0
    kmax = max(k1, k2)
    Cv = view(C, 1:kmax, 1:kmax)
    fill!(Cv, 0)
    for i in 1:kmax
        ni = i <= k1 ? count(blob_a[i]) : 0
        for j in 1:kmax
            if j <= k2
                nj = count(blob_b[j])
                isect = i <= k1 ? count(t -> blob_a[i][t] & blob_b[j][t], 1:N) : 0
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
    blobtransferdistance(blob_a, blob_b)

"""
function blobtransferdistance(
    blob_a::AbstractVector{NTuple{N,Bool}},
    blob_b::AbstractVector{NTuple{N,Bool}},
) where N
    kmax = max(length(blob_a), length(blob_b))
    C = Matrix{Int}(undef, kmax, kmax)
    return blobtransferdistance!(C, blob_a, blob_b)
end