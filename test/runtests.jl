using ADI
using Test
using Statistics
using Random
Random.seed!(8799)

@testset "ADI.jl" begin
    include("decomposition.jl")
end
