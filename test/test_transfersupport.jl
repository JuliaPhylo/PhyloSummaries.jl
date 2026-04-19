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
