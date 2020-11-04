using ADI
using BenchmarkTools
using CSV
using DataFrames
using HCIDatasets: BetaPictoris, HR8799
using InvertedIndices
using PyCall

results = []

# julia benchmarks

julia_algs = (Median(), PCA(20), NMF(20))

for dataset in (BetaPictoris, HR8799)
    @info "Benchmarking - $dataset ADI Reduction"
    cube = dataset[:cube]
    angles = dataset[:pa]
    for (alg, name) in zip(julia_algs, ("median", "pca_20", "nmf_20"))
        (alg isa NMF && dataset === HR8799) && continue
        @info "Julia - $alg"

        time_elapsed = @belapsed $alg($cube, $angles)
        
        @info "Julia - $alg" time=time_elapsed
        push!(results, (framework="ADI.jl", alg=name, N=length(cube), time=time_elapsed))
    end
end

# python benchmarks

vip = pyimport("vip_hci")
py_algs = (
    vip.medsub.medsub,
    (args...; kwargs...) -> vip.pca.pca(args...; ncomp=20,kwargs...),
    (args...; kwargs...) -> vip.nmf.nmf(args...; ncomp=20, kwargs...),
)

for dataset in (BetaPictoris, HR8799)
    @info "Benchmarking - $dataset ADI Reduction"
    cube = dataset[:cube]
    angles = dataset[:pa]
    for (alg, name) in zip(julia_algs, ("median", "pca_20", "nmf_20"))
        (name == "nmf_20" && dataset === HR8799) && continue
        @info "Python - $alg"

        time_elapsed = @belapsed $alg($cube, $angles)
        
        @info "Python - $alg" time=time_elapsed
        push!(results, (framework="VIP", alg=name, N=length(cube), time=time_elapsed))
    end
end

# load results file and update
path = joinpath(@__DIR__, "adi_benchmarks.csv")
df = isfile(path) ? CSV.File(path) |> DataFrame : DataFrame()
out = append!(DataFrame(results), df)
unique!(out, Not(:time)) |> CSV.write(path)
