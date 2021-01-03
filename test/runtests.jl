using ADI
using HCIToolbox
using Test
using Statistics
using LinearAlgebra
using StableRNGs

rng = StableRNG(8799)

@testset "Classic" begin include("classic.jl") end
@testset "PCA" begin include("pca.jl") end
@testset "TPCA" begin include("tpca.jl") end
@testset "NMF" begin include("nmf.jl") end
@testset "GreeDS" begin include("greeds.jl") end

@testset "Metrics" begin 
    include("metrics.jl")
    include("contrast.jl")
end

@testset "SDI" begin include("sdi.jl") end
@testset "Geometries" begin include("geometries.jl") end

