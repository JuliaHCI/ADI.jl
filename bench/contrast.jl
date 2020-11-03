using ADI
using BenchmarkTools
using CSV
using DataFrames
using HCIDatasets: BetaPictoris, HR8799
using HCIToolbox
using PyCall

results = []

# parameters
nbranch = 3
fwhm = 8

# julia benchmarks

alg = PCA(20)
psf = construct(Kernels.Normal(5), (21, 21))

for dataset in (BetaPictoris, HR8799)
    @info "Benchmarking - $dataset Contrast Curve"
    cube = dataset[:cube]
    angles = dataset[:pa]

    time_elapsed = @belapsed contrast_curve($alg, $cube, $angles, $psf; fwhm=fwhm, nbranch=nbranch)

    @info "Julia" time=time_elapsed
    push!(results, (framework="ADI.jl", alg="pca_20", N=length(cube), time=time_elapsed))
end

# python benchmarks

vip = pyimport("vip_hci")

for dataset in (BetaPictoris, HR8799)
    @info "Benchmarking - $dataset Contrast Curve"
    cube = dataset[:cube]
    angles = dataset[:pa]

    time_elapsed = @belapsed vip.metrics.contrast_curve($cube, $angles,  $psf, fwhm=fwhm, pxscale=1e-3, starphot=Metrics.estimate_starphot($cube, fwhm), algo=vip.pca.pca, nbranch=nbranch, verbose=false, plot=false, ncomp=20)

    @info "Python" time=time_elapsed
    push!(results, (framework="VIP", alg="pca_20", N=length(cube), time=time_elapsed))
end

# load results file and update
path = joinpath(@__DIR__, "contrast_benchmarks.csv")
df = isfile(path) ? CSV.File(path) |> DataFrame : DataFrame()
out = append!(DataFrame(results), df)
unique!(out) |> CSV.write(path)
