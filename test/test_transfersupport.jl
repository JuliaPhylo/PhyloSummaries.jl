@testset "blobtransferdistance" begin


β = [(true,false,false,false,false,false,false),
     (false,true,false,false,false,false,false),
     (false,false,true,false,false,false,false),
     (false,false,false,true,false,false,false),
     (false,false,false,false,true,true,true)]

β_star = [(true,false,false,false,false,false,false),
          (false,true,false,false,false,false,false),
          (false,false,true,true,false,false,false),
          (false,false,false,false,true,true,true)]

@test PS.blobtransferdistance(β, β)      == 0
@test PS.blobtransferdistance(β, β_star) == 2
@test PS.blobtransferdistance(β_star, β) == 2

end

@testset "blobtransfer_support, 5 taxa 2 cycles" begin

nfile = joinpath(@__DIR__, "..", "test", "level1_5taxa_2cycles_1blob2bicomp.nwk")
nets = readmultinewick(nfile)
refnet = deepcopy(nets[2])

# 4 of 5 networks have the same blob partition as refnet (dist=0),
# net 3 has blob A,D|C|B|E (dist=2). Total ti_sum = 0+0+2+0+0 = 2
res = @test_logs((:warn, r"non-binary"), (:warn,),
    PS.blobtransfer_support(nets, refnet))
@test length(res.refblobs) == 1
@test res.nnet == 5
@test res.transfer_index_sum == [2.0]

end
