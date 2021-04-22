
cube, angles = BetaPictoris[:cube, :pa]
cube .-= minimum(cube) # rescale so min is > 0
av = AnnulusView(cube; inner=10, outer=20)
mav = MultiAnnulusView(cube, 5)

@testset "Annulus - $(nameof(typeof(alg)))" for alg in [PCA(10), Classic(), NMF(3), GreeDS(3)]
    res = alg(av, angles)
    @test size(res) == (101, 101)
    if !(alg isa NMF)
        res_rdi = alg(av, angles; ref=av)
        @test res ≈ res_rdi
        if alg isa GreeDS
            @test_throws TypeError alg(av, angles; ref=cube)
        else
            @test_throws ErrorException alg(av, angles; ref=cube)
        end
    end

    if alg isa GreeDS
        S = ADI.fit(alg, av; angles=angles) |> reconstruct
    else
        S = ADI.fit(alg, av) |> reconstruct
    end
    X = av()
    @test size(S) == size(X) != size(flatten(cube))
end

# TODO can't test NMF because it's non-deterministic somehow????
@testset "MultiAnnulus - $(nameof(typeof(alg)))" for alg in [PCA(10), Classic()]
    res = alg(mav, angles)
    @test size(res) == (101, 101)
    res_rdi = alg(mav, angles; ref=mav)
    @test res ≈ res_rdi
    @test_throws ErrorException alg(mav, angles; ref=cube)

    des = ADI.fit(alg, mav)
    N = length(mav.indices)
    recons = reconstruct(des)
    @test length(des) == length(recons) == N
    S = reconstruct(alg, mav)
    @test S ≈ inverse(mav, recons)

    # test vector algs
    S2 = reconstruct(fill(alg, N), mav)
    @test S ≈ S2
    S3 = reconstruct(fill(alg, N), mav; ref=mav)
    @test S ≈ S3
    
    # make sure ref is same type
    @test_throws ErrorException reconstruct(fill(alg, N), mav; ref=cube)
end

@testset "av - framewise - $(nameof(typeof(alg)))" for alg in [PCA(10), NMF(2), Classic()]
    fr_alg = Framewise(alg)
    # gets radius automatically
    S1 = reconstruct(fr_alg, av; angles=angles, fwhm=4.7)
    S2 = reconstruct(fr_alg, av; angles=angles, r=15, fwhm=4.7)
    @test S1 ≈ S2
end

@testset "mav - framewise - $(nameof(typeof(alg)))" for alg in [PCA(10), NMF(2), Classic()]
    fr_alg = Framewise(alg)
    # gets fwhm automatically
    S1 = reconstruct(fr_alg, mav; angles=angles)
    S2 = reconstruct(fr_alg, mav; angles=angles, fwhm=5)
    @test S1 ≈ S2

    # test repeated
    k = length(mav.indices)
    algs = Framewise(fill(alg, k))
    S3 = reconstruct(algs, mav; angles=angles)
    @test S1 ≈ S3

    # test delta rot cases
    fr_alg = Framewise(alg, delta_rot=(0.1, 1))
    fr_alg2 = Framewise(alg, delta_rot=range(0.1, 1, length=k))
    S1 = reconstruct(fr_alg, mav; angles=angles)
    S2 = reconstruct(fr_alg2, mav; angles=angles)
    @test S1 ≈ S2
end