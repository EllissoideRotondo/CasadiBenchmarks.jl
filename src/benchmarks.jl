# Imports
using Revise
using LinearAlgebra
using StaticArrays
using FastClosures
using ForwardDiff
using Enzyme
using BenchmarkTools
includet("casadi_min_time.jl")
includet("autodiff_min_time.jl")
includet("modified_equinoctial.jl")

# Test data
x = @SVector rand(14)
thr = 0.1
ceff = 400
μ = 1

# Validate function outputs
f0 = casadi_min_time(x, thr, ceff, μ)
f1 = autodiff_min_time(x, thr, ceff, μ, "ForwardDiff")
f2 = autodiff_min_time(x, thr, ceff, μ, "Enzyme")
@assert f0 ≈ f1
@assert f0 ≈ f2

# Run benchmark
b0 = @benchmark casadi_min_time($x, $thr, $ceff, $μ)
b1 = @benchmark autodiff_min_time($x, $thr, $ceff, $μ, "ForwardDiff")
b2 = @benchmark autodiff_min_time($x, $thr, $ceff, $μ, "Enzyme")
display(b0)
display(b1)
display(b2)