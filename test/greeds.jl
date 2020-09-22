
@testset "Interface" begin
    @test GreeDS().alg == PCA()
    @test GreeDS().threshold == 0
    @test GreeDS(7).alg == PCA(7)
end

@testset "Decomposition" begin
    data = 4 .* randn(rng, 30, 101, 101) .+ 10
    angles = sort!(90randn(rng, 30)) |> normalize_par_angles

    # get sizes correct for ncomps
    for N in [1, 3, 5]
        A, w = @inferred decompose(GreeDS(PCA(ncomps=N)), data, angles)
        @test size(A) == (N, 101 * 101)
        @test size(w) == (N, 30)
    end
    @test_throws ErrorException decompose(GreeDS(PCA(40)), data, angles)
    A, w = decompose(PCA(30; pratio=0.5), data, angles)
    @test size(A, 1) < 30
    @test size(w, 1) < 30

    # default is to use whole cube
    S = reconstruct(GreeDS(PCA()), data, angles)
    @test size(S) == (30, 101, 101)
end

@testset "ADI Trivial" begin
    cube = ones(10, 100, 100) 
    angles = zeros(10)

    reduced_5 = GreeDS(PCA(5))(cube, angles)
    reduced_10 = GreeDS(PCA(10))(cube, angles)

    @test size(reduced_5) == size(reduced_10) == (100, 100)
    @test all(x -> isapprox(x, 0, atol = 1e-9), reduced_5)
    @test all(x -> isapprox(x, 0, atol = 1e-9), reduced_10)
end

@testset "RDI Trivial" begin
    data = randn(rng, 30, 101, 101)
    angles = sort!(90rand(rng, 30)) |> normalize_par_angles

    # not supported for GreeDS alg
    S = data .- reconstruct(GreeDS(PCA(1)), data, angles, zeros(10, 101, 101))
    @test S ≈ data rtol=2e-1
end

@testset "Decomposition - TPCA" begin
    data = 4 .* randn(rng, 30, 101, 101) .+ 10
    angles = sort!(90randn(rng, 30)) |> normalize_par_angles

    # get sizes correct for ncomps
    for N in [1, 3, 5]
        # A, w = @inferred decompose(GreeDS(TPCA(N)), data, angles)
        A, w = decompose(GreeDS(TPCA(N)), data, angles)
        @test size(A) == (N, 101 * 101)
        @test size(w) == (N, 30)
    end
    @test_throws ErrorException decompose(GreeDS(TPCA(40)), data, angles)

    # default is to use whole cube
    S = reconstruct(GreeDS(TPCA()), data, angles)
    @test size(S) == (30, 101, 101)
end

@testset "Decomposition - NMF" begin
    data = 4 .* randn(rng, 10, 101, 101) .+ 10
    angles = sort!(90randn(rng, 10)) |> normalize_par_angles

    # get sizes correct for ncomps
    N = 3
    A, w = decompose(GreeDS(NMF(N)), data, angles)
    @test size(A) == (N, 101 * 101)
    @test size(w) == (N, 10)

    @test_throws ErrorException decompose(GreeDS(NMF(40)), data, angles)

    # default is to use whole cube
    S = reconstruct(GreeDS(NMF()), data, angles)
    @test size(S) == (10, 101, 101)
end
