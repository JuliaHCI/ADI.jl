using ADI
using HCIToolbox
using Test
using Statistics
using LinearAlgebra
using StableRNGs

rng = StableRNG(8799)

@testset "Median" begin include("median.jl") end
@testset "PCA" begin include("pca.jl") end

# @testset "ADI.jl" begin
#     include("pca.jl")
#     include("klip.jl")
# end
