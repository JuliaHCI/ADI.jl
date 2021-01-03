
cube = 4 .* randn(rng, 30, 101, 101) .+ 1000
angles = sort!(90rand(rng, 30))
av = AnnulusView(cube, inner=10, outer=20)
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
    N = length(eachannulus(mav))
    recons = reconstruct(des)
    @test length(des) == length(recons) == N
    S = reconstruct(ALG, mav)
    @test S ≈ inverse(mav, recons)
end