
cube, angles = BetaPictoris[:cube, :pa]

@testset "Interface" begin
    @test GreeDS().kernel == PCA()
    @test GreeDS().threshold == 0
    @test GreeDS(7).kernel == PCA(7)
end

@testset "Decomposition" begin

    # get sizes correct for ncomps
    for N in [1, 3, 5]
        A, w = @inferred ADI.fit(GreeDS(PCA(ncomps=N)), cube; angles=angles)
        @test size(A) == (N, 101 * 101)
        @test size(w) == (size(cube, 1), N)
    end

    # default is to use whole cube
    S = reconstruct(GreeDS(), cube; angles=angles)
    @test size(S) == size(cube)
end

@testset "ADI Trivial" begin
    data = ones(10, 100, 100) 
    angs = zeros(10)

    reduced_5 = GreeDS(PCA(5))(data, angs)
    reduced_10 = GreeDS(PCA(10))(data, angs)

    @test size(reduced_5) == size(reduced_10) == (100, 100)
    @test all(x -> isapprox(x, 0, atol = 1e-9), reduced_5)
    @test all(x -> isapprox(x, 0, atol = 1e-9), reduced_10)
end

@testset "RDI Trivial" begin

    S = subtract(GreeDS(1), cube; angles=angles, ref=zero(cube))
    @test S ≈ cube rtol=2e-1
end

@testset "Decomposition - TPCA" begin
    # get sizes correct for ncomps
    for N in [1, 3, 5]
        A, w = ADI.fit(GreeDS(TPCA(N)), cube; angles=angles)
        @test size(A) == (N, 101 * 101)
        @test size(w) == (size(cube, 1), N)
    end

    # default is to use whole cube
    S = reconstruct(GreeDS(TPCA()), cube; angles=angles)
    @test size(S) == size(cube)
end

# @testset "Decomposition - NMF" begin
#     data = 4 .* randn(rng, 10, 101, 101) .+ 10
#     angles = sort!(90randn(rng, 10)) |> normalize_par_angles

#     # get sizes correct for ncomps
#     N = 3
#     A, w = decompose(GreeDS(NMF(N)), data, angles)
#     @test size(A) == (N, 101 * 101)
#     @test size(w) == (10, N)

#     @test_throws ErrorException decompose(GreeDS(NMF(40)), data, angles)

#     # default is to use whole cube
#     S = reconstruct(GreeDS(NMF()), data, angles)
#     @test size(S) == (10, 101, 101)
# end
