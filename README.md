# SIMDPoly.jl

SIMD optimized polyomial evaluation for real and complex inputs. 

## Usage

The functions for packing and evaluating polynomials are suffixed
with `(N)x(L)[r|c]`:
- `N` is the number of polynomials to evaluate
- `L` is the number of SIMD lanes per polynomial
- `r` or `c` signifies that the polynomial will be evaluated
with a real or complex input, respectively.

For instance evaluating one polynomial with four SIMD lanes:

```julia
coeffs = ntuple(x -> Float64(x), 30) # p(x) = 1 + 2x^1 + ... + 30x^29
packed = packpoly1x4r(coeffs)        # pack coefficients into SIMD vectors
result = evalpoly1x4r(2.0, packed)   # evaluate p(2)
```

## Benchmarks

todo
