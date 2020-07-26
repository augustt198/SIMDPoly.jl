module SIMDPoly

using SIMD

export evalpolydiv,
    pack_coeffs,
    simdrepoly4,
    packrepoly4

include("realpoly.jl")
include("complexpoly.jl")

end
