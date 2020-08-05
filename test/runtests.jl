using ADI
using HCIToolbox
using Test
using Statistics
using LinearAlgebra
using StableRNGs

rng = StableRNG(8799)

@testset "Median" begin include("median.jl") end
@testset "PCA" begin include("pca.jl") end
@testset "Pairet" begin include("pairet.jl") end

