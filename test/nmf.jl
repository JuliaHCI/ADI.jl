
cube, angles = BetaPictoris[:cube, :pa]

@testset "Interface" begin
    @test NMF(5) == NMF(ncomps=5)
    @test isnothing(NMF().ncomps)
end

@testset "Decomposition" begin
    # get sizes correct for ncomps
    for N in [1, 3, 5]
        A, w = ADI.fit(NMF(ncomps=N), cube)
        @test size(A) == (N, 101 * 101)
        @test size(w) == (size(cube, 1), N)
    end

    # default is to use whole cube
    S = reconstruct(NMF(), cube)
    @test size(S) == size(cube)
end

@testset "ADI Trivial" begin
    data = ones(10, 100, 100)
    angs = zeros(10)

    reduced_5 = NMF(5)(data, angs)
    reduced_10 = NMF(10)(data, angs)

    @test size(reduced_5) == size(reduced_10) == (100, 100)
    @test all(x -> isapprox(x, 0, atol = 1e-9), reduced_5)
    @test all(x -> isapprox(x, 0, atol = 1e-9), reduced_10)
end

@testset "RDI Trivial" begin
    S = subtract(NMF(1), cube; ref=zero(cube))
    targ = cube .- minimum(cube)
    @test S ≈ targ rtol=2e-2
end
