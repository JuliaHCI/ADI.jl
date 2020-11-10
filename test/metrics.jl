function test_method(func)
    frame = zeros(101, 101)
    frame[77, 77] = 1.0
    # central fwhm is all NaN
    @test func(frame, (51, 51), 3) === NaN
    # centered on point is infinite, because bkg stddev is 0
    @test func(frame, (77, 77), 3) ≈ Inf
end

@testset "snr" begin
    test_method(snr)
end

@testset "significance" begin
    test_method(significance)
    
    X = randn(rng, 101, 101) .+ 10
    # expect small sample to matter ~3fwhm out
    @test abs(significance(X, (58, 51), 3)) < abs(snr(X, (58, 51), 3))
end

@testset "noise" begin
    data = ones(101, 101)

    @test noise(data, (51, 51), 10) === NaN
end

@testset "detectionmap - $func" for func in (snr, significance, noise)
    data = ones(101, 101)
    _map = detectionmap(func, data, 10)
    # center should be fill value
    @test _map[51, 51] ≈ 0
    # there should be some values for the rest
    @test maximum(_map) > 0
    @test abs(minimum(_map)) ≥ 0
end
