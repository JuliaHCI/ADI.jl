
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
        @test size(A) == (101 * 101, N)
        @test size(w) == (N, size(cube, 3))
    end

    # default is to use whole cube
    S = reconstruct(GreeDS(), cube; angles=angles)
    @test size(S) == size(cube)
end

@testset "ADI Trivial" begin
    data = ones(100, 100, 10)
    angs = zeros(10)

    reduced_5 = GreeDS(PCA(5))(data, angs)
    reduced_10 = GreeDS(PCA(10))(data, angs)

    @test size(reduced_5) == size(reduced_10) == (100, 100)
    @test all(x -> isapprox(x, 0, atol=1e-9), reduced_5)
    @test all(x -> isapprox(x, 0, atol=1e-9), reduced_10)
end

@testset "RDI Trivial" begin

    S = subtract(GreeDS(1), cube; angles=angles, ref=zero(cube))
    @test S ≈ cube rtol=2e-1
end

@testset "Decomposition - NMF" begin
    data = ones(100, 100, 10)
    angs = zeros(10)

    reduced_8 = GreeDS(NMF(8))(data, angs)
    reduced_10 = GreeDS(NMF(10))(data, angs)

    @test size(reduced_8) == size(reduced_10) == (100, 100)
    @test all(x -> isapprox(x, 0, atol=1e-8), reduced_8)
    @test all(x -> isapprox(x, 0, atol=1e-6), reduced_10)
end
