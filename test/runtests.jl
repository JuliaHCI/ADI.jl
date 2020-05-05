using ADI
using HCIToolbox
using Test
using Statistics
using LinearAlgebra
using Random
Random.seed!(8799)

@testset "ADI.jl" begin
    include("pca.jl")
end
