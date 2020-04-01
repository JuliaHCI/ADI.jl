@testset "PCA" begin
    
    @testset "Trivial" begin
        data = ones(30, 512, 512)
        angles = zeros(30)

        res = pca(data, angles)
        @test size(res) == (512, 512)
        @test res ≈ ones(512, 512)
    end

    @testset "RDI Trivial" begin
        data = randn(30, 512, 512)
        angles = 90 .* rand(30)
        
        @test pca(data, angles, zeros(20, 512, 512)) ≈ data

        # interface
        @test @inferred(pca(data, angles)) == @inferred(pca(data, angles, data))
    end
end
