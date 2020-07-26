"""
    evalpoly1x2r

Evaluate one polynomial at `x` given packed coefficients `poly` using
Horner's method with two SIMD lanes.
"""
@generated function evalpoly1x2r(x::T, poly::NTuple{N, Vec{2, T}}) where {N, T<:Real}
    insns = [:(accum = muladd(x2vec, accum, poly[$i])) for i in N-1:-1:1]
    Expr(:block,
        quote
            x2 = x*x
            x2vec = Vec{2, T}((x2, x2))
            accum = poly[N]
        end,
        insns...,
        quote
            accum[1] + accum[2]*x
        end
    )
end

"""
    packpoly1x2r

Pack a tuple of polynomial coefficients into 2-wide SIMD vectors
"""
function packpoly1x2r(p::NTuple{N, T}) where {N, T<:Real}
    parr = collect(p)
    # pad with zeros so length is multiple of 2
    addl = 1 - (N + 1) % 2
    for i = 1:addl
        push!(parr, zero(T))
    end

    m = (N + 1) ÷ 2
    packedarr = Array{Vec{2, T}, 1}(undef, m)
    for i = 1:m
        ofs = (i-1)*2
        packedarr[i] = Vec{2, T}((parr[1+ofs], parr[2+ofs]))
    end
    return tuple(packedarr...)
end

"""
    evalpoly1x4r

Evaluate one polynomial at `x` given packed coefficients `poly` using
Horner's method with four SIMD lanes.
"""
@generated function evalpoly1x4r(x::T, poly::NTuple{N, Vec{4, T}}) where {N, T<:Real}
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
            accum[1] + accum[2]*x + accum[3]*x2 + accum[4]*x2*x
        end
    )
end

"""
    packpoly1x4r

Pack a tuple of polynomial coefficients into 4-wide SIMD vectors
"""
function packpoly1x4r(p::NTuple{N, T}) where {N, T<:Real}
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
