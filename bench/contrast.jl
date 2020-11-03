using ADI
using BenchmarkTools
using CSV
using DataFrames
using HCIDatasets: BetaPictoris, HR8799
using HCIToolbox
using PyCall

results = []

# julia benchmarks

alg = PCA(20)
psf = construct(Kernels.Normal(5), (21, 21))
nbranch = 3

for dataset in (BetaPictoris, HR8799)
    @info "Benchmarking - $dataset Contrast Curve"
    cube = dataset[:cube]
    angles = dataset[:pa]

    time_elapsed = @belapsed contrast_curve($alg, $cube, $angles, $psf; fwhm=5, nbranch=nbranch)

    @info "Julia" time=time_elapsed
    push!(results, (framework="ADI.jl", alg="pca_20", N=length(cube), time=time_elapsed))
end

# python benchmarks

vip = pyimport("vip_hci")
py_algs = (
    vip.medsub.medsub,
    (args...; kwargs...) -> vip.pca.pca(args...; ncomp=20,kwargs...),
    (args...; kwargs...) -> vip.nmf.nmf(args...; ncomp=20, kwargs...),
)

for dataset in (BetaPictoris, HR8799)
    @info "Benchmarking - $dataset Contrast Curve"
    cube = dataset[:cube]
    angles = dataset[:pa]

    time_elapsed = @belapsed vip.metrics.contrast_curve($cube, $angles,  $psf, 5, pxscale=1e-3, starphot=Metrics.starphot($cube, 5), vip.pca.pca, nbranch=nbranch, verbose=false, ncomp=20)

    @info "Python" time=time_elapsed
    push!(results, (framework="VIP", alg="pca_20", N=length(cube), time=time_elapsed))
end

# load results file and update
path = joinpath(@__DIR__, "contrast_benchmarks.csv")
df = isfile(path) ? CSV.File(path) |> DataFrame : DataFrame()
out = append!(DataFrame(results), df)
unique!(out) |> CSV.write(path)
