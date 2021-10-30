using Distances

cube, angles = BetaPictoris[:cube, :pa]

@testset "Interface" begin
    # default args
    @test LOCI() == LOCI(dist_threshold=nothing, metric=Cityblock())
end

@testset "Decomposition" begin
    # get sizes correct for ncomps
    A, w = @inferred ADI.fit(LOCI(), cube)
    @test A ≈ flatten(cube)

    @test size(w) == (size(cube, 1), size(cube, 1))

    # no effect without Framewise
    A2, w2 = ADI.fit(LOCI(dist_threshold=90), cube)
    @test A ≈ A2
    @test w ≈ w2
end

@testset "ADI Trivial" begin
    data = ones(10, 100, 100) 
    angs = zeros(10)

    reduced_full = LOCI()(data, angs)

    @test size(reduced_full) == (100, 100)
    @test all(x -> isapprox(x, 0, atol = 1e-9), reduced_full)

    betapic_reduced_full = LOCI()(cube, angles)
    betapic_reduced_trunc = Framewise(LOCI())(cube, angles)
    betapic_reduced_trunc2 = Framewise(LOCI(dist_threshold=0.90))(cube, angles)
    betapic_reduced_trunc3 = Framewise(LOCI(dist_threshold=0.90, metric=Euclidean()))(cube, angles)
    @test betapic_reduced_full ≉ betapic_reduced_trunc
    @test betapic_reduced_trunc ≉ betapic_reduced_trunc2
    @test betapic_reduced_trunc ≉ betapic_reduced_trunc3
end

@testset "RDI Trivial" begin    
    S = subtract(LOCI(), cube; ref=zero(cube))
    @test S ≈ cube rtol=2e-2
end
