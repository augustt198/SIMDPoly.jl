function evalpolydiv(z::Complex{T}, p1::NTuple{N, T}, p2::NTuple{N, T}) where {N, T}
    coeffs = pack_coeffs(p1, p2)
    return evalpolydiv(z, coeffs)
end

@generated function evalpolydiv(z::Complex, poly::NTuple{N, Vec{4, T}}) where {N, T}
    innerloop = []
    for i = 2:N-2
        insns = quote
            madd = muladd(Avec, RSvec, Bvec)
            Avec = shufflevector(madd, Val(Avec_shuffle_mask))
            Bvec = shufflevector(madd, poly[$(N-i-1)], Val(Bvec_shuffle_mask))
        end
        push!(innerloop, insns)
    end

    Expr(:block,
        quote # setup vectors
            x, y = real(z), imag(z)
            r = 2x
            s = muladd(x, x, y*y)
            # could get rid of one negation instr if LLVM used vfmaddsub
            RSvec = Vec{4, $T}((r, -s, r, -s))
            Avec, Bvec = poly[$N], poly[$N-2]
            madd  = Vec{4, $T}(0)
            Avec_shuffle_mask = (0, 0, 2, 2)
            Bvec_shuffle_mask = (1, 5, 3, 7)
        end,
        innerloop...,
        quote # some common subexprs if we complex divide manually
            # consider special case for leading coefficients being 1 & 0?
            madd = muladd(Avec, RSvec, Bvec) # one left over from inner loop
            a, b, c, d = madd[1], madd[2], madd[3], madd[4]
            denom = muladd(s*c + 2*d*x, c, d*d)
            numre = muladd(s*c, a, x*(a*d + b*c) + b*d)
            numim = y*(a*d - b*c)
            return Complex(numre/denom, numim/denom)
        end
    )
end

# pack two polynomial coefficients into SIMD vectors
function pack_coeffs(p1::NTuple{N, T}, p2::NTuple{N, T}) where {N, T}
    arr = Array{Vec{4, T}, 1}(undef, 0)
    for i = 1:N
        # the index N-2 vector of coefficients is modified so that it can also be
        # used as the initial B vector in `evalpolydiv`, to reduce setup time
        idxoff = (i == N-2) ? 1 : 0
        vec = Vec{4, T}((p1[i+idxoff], p1[i], p2[i+idxoff], p2[i]))
        push!(arr, vec)
    end
    return tuple(arr...)
end
