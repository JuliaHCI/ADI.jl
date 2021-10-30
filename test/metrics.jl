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

@testset "stim map" begin
    S = randn(rng, 100, 512, 512)
    angles = sort!(90rand(rng, 100))

    sm = stimmap(S, angles)
    @test size(sm) == (size(S, 2), size(S, 3))
    @test all(isfinite, sm)
    @test eltype(sm) == eltype(S)

    smt = stim_threshold(sm, S, angles)
    smt2 = stim_threshold(S, angles)
    @test smt ≈ smt2
    @test smt ≈ 0.5 rtol=1e-2

end

@testset "slim map" begin
    cube, angles = BetaPictoris[:cube, :pa]

    resid_cubes = subtract.(PCA.(5:10), Ref(cube))
    stim_maps = stimmap.(resid_cubes, Ref(angles))
    # 1
    stim_av, mask = slimmap(stim_maps; N=20)
    slim1 = stim_av .* mask
    #2
    stim_av2, mask2 = slimmap(resid_cubes, angles; N=20)
    slim2 = stim_av2 .* mask2

    @test all(slim1 .≈ slim2)

    # make sure beta pic b was found
    @test slim1[62, 62] > 0.5
end
