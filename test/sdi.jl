
@testset "SDI Interface" begin
    @test DoubleSDI(PCA(15)) == DoubleSDI(PCA(15), PCA(15))
    @test SliceSDI(PCA(15)) == SliceSDI(PCA(15), PCA(15))
    cube = ones(100, 100, 5, 10)
    angles = zeros(10)
    scales = collect(range(1.2, 1, length=5))
    @test PCA(15)(cube, angles, scales) == SingleSDI(PCA(15))(cube, angles, scales)
end

@testset "Trivial SDI - $(typeof(ALG))" for ALG in (PCA(5), GreeDS(5), Classic(), LOCI())
    cube = ones(100, 100, 5, 10)
    angles = zeros(10)
    scales = collect(range(1.2, 1, length=5))

    reduced= ALG(cube, angles, scales)

    @test size(reduced) == (100, 100)
    @test all(x -> isapprox(x, 0, atol = 1e-9), reduced)
end

@testset "Trivial Double SDI - $(typeof(ALG1)), $(typeof(ALG2))" for ALG1 in (PCA(2), Classic(), GreeDS(2), LOCI()),
                                                 ALG2 in (PCA(5), GreeDS(5), Classic(), LOCI())
    cube = ones(100, 100, 5, 10)
    angles = zeros(10)
    scales = collect(range(1.2, 1, length=5))

    reduced = DoubleSDI(ALG1, ALG2)(cube, angles, scales)

    @test size(reduced) == (100, 100)
    @test all(x -> isapprox(x, 0, atol = 1e-9), reduced)
end

@testset "Trivial Slice SDI - $(typeof(ALG1)), $(typeof(ALG2))" for ALG1 in (PCA(2), Classic(), GreeDS(2), LOCI()),
                                                 ALG2 in (PCA(5), GreeDS(5), Classic(), LOCI())
    cube = ones(100, 100, 5, 10)
    angles = zeros(10)
    scales = collect(range(1.2, 1, length=5))

    reduced = SliceSDI(ALG1, ALG2)(cube, angles, scales)

    @test size(reduced) == (100, 100)
    @test all(x -> isapprox(x, 0, atol = 1e-9), reduced)
end