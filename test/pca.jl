@testset "PCA" begin
    
    @testset "Interface" begin
        data = 4 .* randn(30, 512, 512) .+ 10
        angles = sort!(90randn(30)) |> normalize_par_angles

        # get sizes correct for ncomps
        for N in [1, 3, 5]
            res = pca(data, angles; ncomps=N)
            @test size(res.A) == (N, 512 * 512)
            @test size(res.w) == (N, 30)
        end
        d = pca(data, angles, ncomps=40)
        @test size(d.A) == (30, 512 * 512)
        @test size(d.w) == (30, 30)
        d = pca(data, angles, ncomps=40, pratio=100)
        @test size(d.A) == (30, 512 * 512)
        @test size(d.w) == (30, 30)
        d = pca(data, angles, ncomps=10, pratio=100)
        @test size(d.A) == (10, 512 * 512)
        @test size(d.w) == (10, 30)
        d = pca(data, angles, ncomps=30, pratio=0.5)
        @test size(d.A, 1) < 30
        @test size(d.w, 1) < 30

        @test size(d.S) == (30, 512, 512)
        # angles are the same
        @test d.angles == angles
        @test @inferred(reduce(d)) ≈ @inferred(reduce(d, data)) ≈ @inferred(reduce(d, data, angles))
    end

    @testset "RDI Trivial" begin
        data = randn(30, 512, 512)
        angles = sort!(90randn(30)) |> normalize_par_angles
        
        design = pca(data, zeros(size(data)), angles; ncomps=1)
        @test design.S ≈ data rtol=1e-2
    end
end
