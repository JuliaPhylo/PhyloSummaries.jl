"""
    referencesupport(networks, referencenet; minimumblobdegree=4)

Calculate the support for features (blob partitions, bipartitions, circular
orders, hybrid clades) present in a reference network `referencenet`, based on
their frequency in a sample of `networks`.

"""
function referencesupport(
    networks::AbstractVector{PN.HybridNetwork},
    referencenet::PN.HybridNetwork;
    minimumblobdegree::Int=4,
)

end

