
@testset "evaluation" begin
    cube = randn(rng, 10, 101, 101) .+ 100
    angles = 90 .* rand(rng, 10)
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
