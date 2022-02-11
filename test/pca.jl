
cube, angles = BetaPictoris[:cube, :pa]

@testset "Interface" begin
    @test PCA(5) == PCA(ncomps=5)
    @test isnothing(PCA().ncomps)
end

@testset "Decomposition" begin
    # get sizes correct for ncomps
    for N in [1, 3, 5]
        A, w = @inferred ADI.fit(PCA(ncomps=N), cube)
        @test size(A) == (101 * 101, N)
        @test size(w) == (N, size(cube, 1))
    end
    A, w = ADI.fit(PCA(:pratio; pratio=0.5), cube)
    @test size(A, 1) < size(cube, 3)
    @test size(w, 2) < size(cube, 3)

    A, w = ADI.fit(PCA(:noise), cube)
    @test size(A, 1) ≤ size(cube, 3)
    @test size(w, 2) ≤ size(cube, 3)

    # default is to use whole cube
    S = reconstruct(PCA(), cube)
    @test size(S) == size(cube)
end

@testset "ADI Trivial" begin
    data = ones(100, 100, 10)
    angs = zeros(10)

    reduced_5 = PCA(5)(data, angs)
    reduced_10 = PCA(10)(data, angs)

    @test size(reduced_5) == size(reduced_10) == (100, 100)
    @test all(x -> isapprox(x, 0, atol = 1e-9), reduced_5)
    @test all(x -> isapprox(x, 0, atol = 1e-9), reduced_10)
end

@testset "RDI Trivial" begin
    S = subtract(PCA(1), cube; ref=zero(cube))
    @test S ≈ cube rtol=2e-2
end
