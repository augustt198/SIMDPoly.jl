using SIMDPoly
using Test

@testset "complex poly" begin
    @test evalpolydiv(2.0+0.0im, (1.,2.,3.,4.,5.), (2.,4.,6.,8.,10.)) == evalpoly(2.0+0.0im, (1.,2.,3.,4.,5.))/evalpoly(2.0+0.0im+0.0im+0.0im, (2.,4.,6.,8.,10.))
end
