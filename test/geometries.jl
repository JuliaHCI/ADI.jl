
cube, angles = BetaPictoris[:cube, :pa]
av = AnnulusView(cube; inner=10, outer=20)
mav = MultiAnnulusView(cube, 5)

@testset "Annulus - $ALG" for ALG in [PCA(10), Classic(), NMF(3), GreeDS(3)]
    res = ALG(av, angles)
    @test size(res) == (101, 101)
    if !(ALG isa NMF)
        res_rdi = ALG(av, angles; ref=av)
        @test res ≈ res_rdi
        if ALG isa GreeDS
            @test_throws TypeError ALG(av, angles; ref=cube)
        else
            @test_throws ErrorException ALG(av, angles; ref=cube)
        end
    end

    if ALG isa GreeDS
        S = ADI.fit(ALG, av; angles=angles) |> reconstruct
    else
        S = ADI.fit(ALG, av) |> reconstruct
    end
    X = av()
    @test size(S) == size(X) != size(flatten(cube))
end

@testset "MultiAnnulus - $ALG" for ALG in [PCA(10), Classic(), NMF(3)]
    res = ALG(mav, angles)
    @test size(res) == (101, 101)
    if !(ALG isa NMF)
        res_rdi = ALG(mav, angles; ref=mav)
        @test res ≈ res_rdi
        @test_throws ErrorException ALG(mav, angles; ref=cube)
    end

    des = ADI.fit(ALG, mav)
    N = length(mav.indices)
    recons = reconstruct(des)
    @test length(des) == length(recons) == N
    S = reconstruct(ALG, mav)
    @test S ≈ inverse(mav, recons)

    # test vector algs
    S2 = reconstruct(fill(ALG, N), mav)
    @test S ≈ S2
    if !(ALG isa NMF)
        S3 = reconstruct(fill(ALG, N), mav; ref=mav)
        @test S ≈ S3
    end
    
    # make sure ref is same type
    @test_throws ErrorException reconstruct(fill(ALG, N), mav; ref=cube)
end

@testset "av - framewise - $alg" for alg in [PCA(10), Classic()]
    fr_alg = Framewise(alg)
    # gets radius automatically
    S1 = reconstruct(fr_alg, av; angles=angles, fwhm=4.7)
    S2 = reconstruct(fr_alg, av; angles=angles, r=15, fwhm=4.7)
    @test S1 ≈ S2
end

@testset "mav - framewise - $alg" for alg in [PCA(10), Classic()]
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