


tree1 = readnewick("((A,B),(C,D));")
tree2 = readnewick("((C,D),(B,A));") # ≈ tree1
tree3 = readnewick("((A,C),(B,D));")
tree4 = readnewick("(D,(C,(A,B)));") # ≈ trees 1,2 if unrooted
tree5 = readnewick("(D, C,(A,B)) ;") # ≈ tree4 but root suppressed
taxa = ["A","B","C","D"]

@testset "count_bipartitions!" begin
# unrooted
counts = Dictionary{BitVector,Int}()
PhyloSummaries.count_bipartitions!(counts, tree1, taxa, false)
PhyloSummaries.count_bipartitions!(counts, tree3, taxa, false)
expected = Dict(
    BitVector([1, 1, 0, 0]) => 1,
    BitVector([1, 0, 1, 0]) => 1,
)
@test length(counts) == length(expected)
for bp in keys(expected)
    @test counts[bp] == expected[bp]
end
# rooted
empty!(counts)
PhyloSummaries.count_bipartitions!(counts, tree1, taxa, true)
PhyloSummaries.count_bipartitions!(counts, tree3, taxa, true)
expected = Dict(
    BitVector([1, 1, 0, 0]) => 1,
    BitVector([0, 0, 1, 1]) => 1,
    BitVector([1, 0, 1, 0]) => 1,
    BitVector([0, 1, 0, 1]) => 1,
)
@test length(counts) == length(expected)
for bp in keys(expected)
    @test counts[bp] == expected[bp]
end
end

@testset "consensustree" begin

@test_throws ArgumentError consensustree(PN.HybridNetwork[])
@test_throws ArgumentError consensustree([readnewick("((A,#H1),(B)#H1);")])
@test_throws "not share the same taxon set" consensustree([tree1, readnewick("((A,B),C);")])
@test_throws "not in taxon list" consensustree([tree1, readnewick("((A,B),C,E);")])

# single-tree input
@test writenewick(consensustree([tree1])) == "(C,D,(A,B));"
@test writenewick(consensustree([tree1]; rooted=true)) == "((A,B),(C,D));"

# 4 taxa, 5 trees, missing edge lengths
trees = [tree2,tree1,tree3,tree4,tree5]
con = consensustree(trees; rooted=true, proportion=0.8)
writenewick(con,round=true) == "(A,B,C,D);"
con = consensustree(trees; proportion=0.7)
writenewick(con,round=true) == "(C,D,(B,A):0.8);"
con = consensustree(trees; rooted=true) # greedy
@test writenewick(con,round=true) == "((D,C):0.4,(B,A):0.8);"
# fixit: currently ((B,A):0.8,(D,(C):0.2):0.4);
# problem: extra degree-2 node to C
@test [n.fvalue for n in con.node if !n.leaf] == [-1,.4,.8]
@test [e.y for e in con.edge if !isexternal(e)] == [.4,.8]

tfile = joinpath(@__DIR__,"..","test","raxmltrees.tre")
# tfile = joinpath(dirname(pathof(PhyloSummaries)), "..","test","raxmltrees.tre")
trees = readmultinewick(tfile)
con = consensustree(trees)
@test writenewick(con,round=true) == "(E,O,((A,B):0.833,(C,D):1.0):0.533);" 
# fixit: problem with 3 extra nodes...

end # of sub-testset

#=
@testset "majority-rule consensus (paper example)" begin
    trees = [tree3, tree2, tree1]

    # Expected majority-rule consensus tree
    # (A,B) appears in 2/3 trees → included
    # (C,D) appears in 2/3 trees → included
    # → consensus should be ((A,B),(C,D));

    consensus = consensustree(trees; rooted=true, proportion=0.5)
    @test writenewick(consensus, round = true) == "((D,C):0.667,(B,A):0.667);"
end
=#
@testset "majority-rule consensus 2" begin
    trees = [tree2, tree1]

    # Expected majority-rule consensus tree
    # (A,B) appears in 2/3 trees → included
    # (C,D) appears in 2/3 trees → included
    # → consensus should be ((A,B),(C,D));

    consensus = consensustree(trees; rooted=true, proportion=0.5)
    @test writenewick(consensus) == "((D,C):1.0,(B,A):1.0);"
end


@testset "consensustree (rooted majority, 5 trees, 5 taxa)" begin
        
    nwk = [
        "(((A,B),(C,D)),E);",
        "((E,(A,B)),(C,D));",
        "((A,B),(C,(D,E)));",
        "(((A,B),E),(C,D));",
        "((A,B),((C,D),E));",
    ]
    trees = readnewick.(nwk)

    con = consensustree(trees; rooted=true, proportion=0.5)


    @test writenewick(con) =="(E,(B,A):1.0,(D,C):0.8);"
end



