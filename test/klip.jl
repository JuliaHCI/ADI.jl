@testset "klip - $ALG" for ALG in [PCA, NMF]
    cube = ones(10, 100, 100) 
    angles = zeros(10)

    reduced_5 = reduce(ALG, cube, angles, 5)
    reduced_10 = reduce(ALG, cube, angles, 10)

    @test size(reduced_5) == size(reduced_10) == (100, 100)
    @test all(isapprox.(reduced_5, 0, atol = 1e-9))
    @test all(isapprox.(reduced_10, 0, atol = 1e-9))
end

# @testset "Interface - $ALG" for ALG in [pca, NMF, Median, Mean]
#     cube = ones(10, 100, 100) 
#     angles = 50rand(10)

#     d = design(ALG, cube)
#     @test d.A * d.w ≈ d.S

#     @test_throws MethodError design(ALG, ones(100, 100))
#     @test_throws MethodError design(ALG, ones(10, 50, 100, 100))
# end

# @testset "Trivial Stats Estimators" begin
#     cube = rand(10, 512, 512)
#     d = design(Mean, cube)
#     @test d.A ≈ d.S ≈ reshape(mean(cube, dims = 1), 1, 262144)

#     d = design(Median, cube)
#     @test d.A ≈ d.S ≈ reshape(median(cube, dims = 1), 1, 262144)
# end