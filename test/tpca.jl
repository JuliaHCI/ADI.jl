@testset "Interface" begin
    @test TPCA(5) == TPCA(ncomps=5)
    @test isnothing(TPCA().ncomps)
end

@testset "Decomposition" begin
    data = 4 .* randn(rng, 30, 101, 101) .+ 10
    angles = sort!(90randn(rng, 30)) |> normalize_par_angles

    # get sizes correct for ncomps
    for N in [1, 3, 5]
        # A, w = @inferred decompose(TPCA(N), data, angles) TODO type stability
        A, w = decompose(TPCA(N), data, angles)
        @test size(A) == (N, 101 * 101)
        @test size(w) == (N, 30)
    end
    @test_throws ErrorException decompose(TPCA(40), data, angles)

    # default is to use whole cube
    S = reconstruct(TPCA(), data, angles)
    @test size(S) == (30, 101, 101)
end

# @testset "ADI Trivial" begin
#     cube = ones(10, 100, 100) 
#     angles = zeros(10)

#     reduced_5 = TPCA(5)(cube, angles)
#     reduced_10 = TPCA(10)(cube, angles)

#     @test size(reduced_5) == size(reduced_10) == (100, 100)
#     @test all(x -> isapprox(x, 0, atol = 1e-9), reduced_5)
#     @test all(x -> isapprox(x, 0, atol = 1e-9), reduced_10)
# end

# @testset "RDI Trivial" begin
#     data = randn(rng, 30, 101, 101)
#     angles = sort!(90rand(rng, 30)) |> normalize_par_angles
    
#     S = data .- reconstruct(TPCA(2), data, angles, zeros(10, 101, 101))
#     @test S ≈ data rtol=2e-2
# end
