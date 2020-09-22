using ADI
using HCIToolbox
using Test
using Statistics
using LinearAlgebra
using StableRNGs

rng = StableRNG(8799)

@testset "Median" begin include("median.jl") end
@testset "PCA" begin include("pca.jl") end
@testset "TPCA" begin include("tpca.jl") end
@testset "NMF" begin include("nmf.jl") end
@testset "GreeDS" begin include("greeds.jl") end

