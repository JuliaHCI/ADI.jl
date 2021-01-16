using ADI
using HCIDatasets: BetaPictoris
using LinearAlgebra
using PSFModels
using StableRNGs
using Statistics
using Test

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
@testset "Framewise" begin include("framewise.jl") end
@testset "Geometries" begin include("geometries.jl") end
