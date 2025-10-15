
@testset "consensus trees" begin

tree1 = readnewick("((A,B),(C,D));")
tree2 = readnewick("((A,C),(B,D));")
taxa = ["A","B","C","D"]

@testset "count_bipartitions! unrooted" begin
counts = Dict{NTuple{4,Int}, Int}()
PhyloSummaries.count_bipartitions!(counts, tree1, taxa, false)
PhyloSummaries.count_bipartitions!(counts, tree2, taxa, false)
expected = Dict(
    (1, 1, 0, 0) => 1,
    (1, 0, 1, 0) => 1,
)
@test counts == expected
end


@testset "count_bipartitions! rooted" begin
counts = Dict{NTuple{4,Int}, Int}()
PhyloSummaries.count_bipartitions!(counts, tree1, taxa, true)
PhyloSummaries.count_bipartitions!(counts, tree2, taxa, true)
expected = Dict(
    (1, 1, 0, 0) => 1,
    (0, 0, 1, 1) => 1,
    (1, 0, 1, 0) => 1,
    (0, 1, 0, 1) => 1,
)
@test counts == expected
end

end
