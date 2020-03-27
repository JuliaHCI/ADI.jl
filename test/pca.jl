@testset "PCA" begin
    
    @testset "Trivial" begin
        data = ones(30, 512, 512)
        angles = zeros(30)

        res = pca(data, angles)
        @test size(res) == (512, 512)
        @test res ≈ ones(512, 512)
    end
end
