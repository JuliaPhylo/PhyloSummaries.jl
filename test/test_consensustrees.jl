
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
PhyloSummaries.consensustree([tree1, tree2];)
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

@testset "tuple_from_clustervector" begin
    v1 = [0, 0, 1, 1]
    v2 = [1, 1, 0, 0]
    v3 = [1, 1, 1, 1]  # trivial
    v4 = [0, 0, 0, 0]  # trivial

    @test PhyloSummaries.tuple_from_clustervector(v1, false) == BitVector([0,0,1,1])
    @test PhyloSummaries.tuple_from_clustervector(v2, false) == BitVector([1,1,0,0])
    @test isnothing(PhyloSummaries.tuple_from_clustervector(v3, false))
    @test isnothing(PhyloSummaries.tuple_from_clustervector(v4, false))
end

@testset "consensus_bipartition filtering" begin
    taxa = ["A","B","C","D"]
    splitcounts = Dict(
        BitVector([1,1,0,0]) => 2,
        BitVector([1,0,1,0]) => 1,
        BitVector([0,1,0,1]) => 3
    )
    result = PhyloSummaries.consensus_bipartition(splitcounts, 0.5, 4, taxa)
    # expect bipartitions with freq >= ceil(0.5*4)=2 → first and third included
    @test any(x -> x == BitVector([1,1,0,0]), result)
    @test any(x -> x == BitVector([0,1,0,1]), result)
    @test !any(x -> x == BitVector([1,0,1,0]), result)
end

@testset "create_tree_from_bipartition_set basic topology" begin
    taxa = ["A","B","C","D"]
    biparts = [BitVector([1,1,0,0])]
    tree = PhyloSummaries.create_tree_from_bipartition_set(taxa, biparts)
    @test occursin("A", writeTopology(tree))
    @test occursin("B", writeTopology(tree))
    @test occursin("C", writeTopology(tree))
    @test occursin("D", writeTopology(tree))
end

@testset "consensustree single-tree copy" begin
    t = readnewick("((A,B),(C,D));")
    result = PhyloSummaries.consensustree([t])
    @test writeTopology(result) == writeTopology(t)
end

@testset "consensustree majority rule" begin
    t1 = readnewick("((A,B),(C,D));")
    t2 = readnewick("((A,C),(B,D));")
    t3 = readnewick("((A,B),(C,D));")  # reinforce first split

    consensus = PhyloSummaries.consensustree([t1, t2, t3])
    newick = writeTopology(consensus)
    # majority split (A,B)|(C,D) should dominate
    @test occursin("(A,B)", newick)
    @test occursin("(C,D)", newick)
end

@testset "consensustree argument errors" begin
    using PhyloNetworks
    t = readnewick("(A,B);")
    bad = deepcopy(t)
    bad.numhybrids = 1
    @test_throws ArgumentError PhyloSummaries.consensustree([])
    @test_throws ArgumentError PhyloSummaries.consensustree([bad])
end

end
