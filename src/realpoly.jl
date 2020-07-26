"""
    simdrepoly4

Evaluate a polynomial at `x` given packed coefficients `poly` using
Horner's method with four SIMD lanes.
"""
@generated function simdrepoly4(x::T, c::T, poly::NTuple{N, Vec{4, T}}) where {N, T<:Real}
    insns = [:(accum = muladd(x4vec, accum, poly[$i])) for i in N-1:-1:1]
    Expr(:block,
        quote
            x2 = x*x
            x4 = x2*x2
            x4vec = Vec{4, T}((x4, x4, x4, x4))
            accum = poly[N]
        end,
        insns...,
        quote
            # modify x4vec to form (x, x², x³, x⁴) vector
            setindex(x4vec, x, 1)
            setindex(x4vec, x2, 2)
            setindex(x4vec, x2*x, 3)
            c + sum(x4vec * accum)
        end
    )
end

"""
    packrepoly4

Pack a tuple of polynomial coefficients into 4-wide SIMD vectors
"""
function packrepoly4(p::NTuple{N, T}) where {N, T<:Real}
    parr = collect(p)
    # pad with zeros so length is multiple of 4
    addl = 3 - (N + 3) % 4
    for i = 1:addl
        push!(parr, zero(T))
    end

    m = (N + 3) ÷ 4
    packedarr = Array{Vec{4, T}, 1}(undef, m)
    for i = 1:m
        ofs = (i-1)*4
        packedarr[i] = Vec{4, T}((parr[1+ofs], parr[2+ofs], parr[3+ofs], parr[4+ofs]))
    end
    return tuple(packedarr...)
end
