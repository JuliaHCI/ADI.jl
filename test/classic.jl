
cube, angles = BetaPictoris[:cube, :pa]

@testset "interface" begin
    X = flatten(cube)
    med = median(X, dims=3)

    d = ADI.fit(Classic(), cube)
    @test ADI.design(d) ≈ med

    t = reconstruct(d)
    @test size(t) == size(X)
    @test X .- t ≈ X .- med

end

@testset "evaluation" begin
    res = Classic()(cube, angles)
    res2 = Classic()(cube, angles, ref=cube)
    @test res == res2
    @test median(res) ≈ 0 atol = 1 / sqrt(length(res)) # within 1σ

    # using different statistic
    res = Classic(mean)(cube, angles)
    res2 = Classic(mean)(cube, angles, ref=cube)
    @test res == res2
    @test median(res) ≈ 0 atol = 1 / sqrt(length(res)) # within 1σ
end
