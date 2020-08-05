
@testset "Interface" begin
    @test Pairet().alg == PCA()
    @test Pairet().threshold == 0

    @test Pairet(PCA(5)).alg == PCA(5)
end

@testset "Decomposition" begin
    data = 4 .* randn(rng, 30, 512, 512) .+ 10
    angles = sort!(90randn(rng, 30)) |> normalize_par_angles

    # get sizes correct for ncomps
    for N in [1, 3, 5]
        A, w = @inferred decompose(Pairet(PCA(ncomps=N)), data, angles)
        @test size(A) == (N, 512 * 512)
        @test size(w) == (N, 30)
    end
    @test_throws ErrorException decompose(Pairet(PCA(40)), data, angles)
    A, w = decompose(PCA(30; pratio=0.5), data, angles)
    @test size(A, 1) < 30
    @test size(w, 1) < 30

    # default is to use whole cube
    S = reconstruct(Pairet(), data, angles)
    @test size(S) == (30, 512, 512)
end

@testset "ADI Trivial" begin
    cube = ones(10, 100, 100) 
    angles = zeros(10)

    reduced_5 = Pairet(PCA(5))(cube, angles)
    reduced_10 = Pairet(PCA(10))(cube, angles)

    @test size(reduced_5) == size(reduced_10) == (100, 100)
    @test all(x -> isapprox(x, 0, atol = 1e-9), reduced_5)
    @test all(x -> isapprox(x, 0, atol = 1e-9), reduced_10)
end

@testset "RDI Trivial" begin
    data = randn(rng, 30, 512, 512)
    angles = sort!(90randn(rng, 30)) |> normalize_par_angles

    # not supported for Pairet alg
    @test_throws MethodError reconstruct(Pairet(PCA(1)), data, angles, zeros(10, 512, 512))
end
