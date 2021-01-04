cube = 4 .* randn(rng, 30, 101, 101) .+ 100
angles = sort!(90rand(30))
psf = ones(11, 11)

@testset "throughput - $alg" for alg in [Classic(), PCA(10), GreeDS(3)]
    tt, meta = throughput(alg, cube, angles, psf; fwhm=4)

    @test length(tt) == 10
    @test length(meta.distance) == length(meta.noise) == 10
    @test size(meta.fake_comps) == size(cube)[2:3]
    @test keys(meta) == (:distance, :fake_comps, :noise)

    t1 = throughput(alg, cube, angles, psf, (51, 61); fwhm=4)
    t2 = throughput(alg, cube, angles, psf, Polar(10, deg2rad(90)); fwhm=4)
    @test t1 ≈ t2
end

@testset "contrast - $alg" for alg in [Classic(), PCA(10), GreeDS(3)]
    cc = contrast_curve(alg, cube, angles, psf; fwhm=4)

    @test keys(cc) == (:distance, :throughput, :contrast, :contrast_corr, :noise)
    ls = map(length, values(cc))
    @test all(ls .== 38)


    reduced_empty = alg(cube, angles)
    sp = Metrics.estimate_starphot(cube, 4)
    cc_raw = contrast_curve(alg, cube, angles, psf; fwhm=4, subsample=false, starphot=sp)
    cc_sub = Metrics.subsample_contrast(reduced_empty, cc_raw.distance, cc_raw.throughput; fwhm=4, starphot=sp)
    
    for (v1, v2) in zip(values(cc), values(cc_sub))
        v1_ = filter(isfinite, v1)
        v2_ = filter(isfinite, v2)
        @test v1_ ≈ v2_
    end
end

@testset "throughput - AnnulusView - $alg" for alg in [Classic(), PCA(10), GreeDS(3)]
    av = AnnulusView(cube; inner=8)
    tt, meta = throughput(alg, av, angles, psf; fwhm=4)

    @test length(tt) == 9
    @test length(meta.distance) == length(meta.noise) == 9
    @test size(meta.fake_comps) == size(cube)[2:3]
    @test keys(meta) == (:distance, :fake_comps, :noise)

    t1 = throughput(alg, cube, angles, psf, (51, 61); fwhm=4)
    t2 = throughput(alg, cube, angles, psf, Polar(10, deg2rad(90)); fwhm=4)
    @test t1 ≈ t2
end

@testset "throughput - MultiAnnulusView - $alg" for alg in [Classic(), PCA(10)]
    av = MultiAnnulusView(cube, 4; inner=8)
    tt, meta = throughput(alg, av, angles, psf; fwhm=4)

    @test length(tt) == 9
    @test length(meta.distance) == length(meta.noise) == 9
    @test size(meta.fake_comps) == size(cube)[2:3]
    @test keys(meta) == (:distance, :fake_comps, :noise)

    t1 = throughput(alg, cube, angles, psf, (51, 61); fwhm=4)
    t2 = throughput(alg, cube, angles, psf, Polar(10, deg2rad(90)); fwhm=4)
    @test t1 ≈ t2
end
