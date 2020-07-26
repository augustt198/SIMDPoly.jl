using SIMDPoly
using BenchmarkTools
using PyPlot

NSAMPLES = 10_000

degrees = 3:40
simd_times, ref_times = [], []

for n in degrees
    poly1 = ntuple(x -> Float64(x), n)
    poly2 = 2 .* poly1
    packed = packpoly2x2c(poly1, poly2)
    our = @benchmark evalpolydiv2x2c(z, $packed) samples=NSAMPLES setup=(z=rand()+rand()*im)
    ref = @benchmark (evalpoly(z, $poly1) / evalpoly(z, $poly2)) samples=NSAMPLES setup=(z=rand()+rand()*im)
    push!(simd_times, minimum(our).time)
    push!(ref_times, minimum(ref).time)
    #println(minimum(our).time, " / ", minimum(ref).time)
end

ioff()
plot(degrees, ref_times, linewidth=3, label="reference (Base.Math.evalpoly)")
plot(degrees, simd_times, linewidth=3, label="SIMD (4 x f64 lanes)")
legend()
ylim(ymin=0)
grid(axis="y", alpha=0.5)
xlabel("N")
ylabel("time (ns)")
title("degree N polynomial divison, complex input")
savefig("bm_plot_complex.svg")
