
@testset "evaluation" begin
    cube = randn(rng, 10, 101, 101) .+ 100
    angles = 90 .* rand(rng, 10)
    res = Median()(cube, angles)
    res2 = Median()(cube, angles, cube)
    @test res == res2
    @test median(res) ≈ 0 atol = 1 / sqrt(length(res)) # within 1σ
end
