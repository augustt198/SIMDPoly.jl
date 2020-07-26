module SIMDPoly

using SIMD

export evalpolydiv2x2c,
    packpoly2x2c,
    evalpoly1x2r,
    evalpoly1x4r,
    packpoly1x2r,
    packpoly1x4r

include("realpoly.jl")
include("complexpoly.jl")

end
