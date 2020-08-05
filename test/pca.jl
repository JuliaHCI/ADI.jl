@testset "PCA" begin
    
    @testset "Interface" begin
        @test PCA(5) == PCA(ncomps=5)
        @test isnothing(PCA().ncomps)
        @test PCA(5).pratio == 1
    end

    @testset "Decomposition" begin
        data = 4 .* randn(rng, 30, 512, 512) .+ 10
        angles = sort!(90randn(rng, 30)) |> normalize_par_angles

        # get sizes correct for ncomps
        for N in [1, 3, 5]
            A, w = @inferred decompose(PCA(ncomps=N), data, angles)
            @test size(A) == (N, 512 * 512)
            @test size(w) == (N, 30)
        end
        @test_throws ErrorException decompose(PCA(40), data, angles)
        A, w = decompose(PCA(30; pratio=0.5), data, angles)
        @test size(A, 1) < 30
        @test size(w, 1) < 30


        S = reconstruct(PCA(), data, angles)
        @test size(S) == (30, 512, 512)
    end

    @testset "RDI Trivial" begin
        data = randn(rng, 30, 512, 512)
        angles = sort!(90randn(rng, 30)) |> normalize_par_angles
        
        A, w = decompose(PCA(1), data, angles, zeros(10, 512, 512))
        @test_broken all(x -> isapprox(x, 0, atol=1e-5), A) # TODO
    end
end
