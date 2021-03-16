using ADI: find_angles, compute_pa_thresh

@testset "helpers" begin
    angles = 10:10:90
    fwhm = 4.7
    r = 15
    idx = 6

    full_inds = axes(angles, 1)
    @test find_angles(angles, idx, 3) == setdiff(full_inds, idx)
    @test find_angles(angles, idx, 3, limit=4) == [4, 5, 7, 8]

    @test compute_pa_thresh(angles, r, fwhm, 1.5) ≈ 26.449102384427057
    @test compute_pa_thresh(angles, r, fwhm, 5) ≈ 36
end

cube, angles = BetaPictoris[:cube, :pa]

@testset "framewise - $alg" for alg in [PCA(10), NMF(3), Classic()]
    fr_alg = Framewise(alg)
    S = reconstruct(fr_alg, cube; angles=angles, fwhm=4.7, r=15)
    @test size(S) == size(cube)

    @test fr_alg(cube, angles; r=15, fwhm=4.7) ≈ collapse!(cube .- S, angles)
end
