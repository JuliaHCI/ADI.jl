using ADI
using BenchmarkTools
using CSV
using DataFrames
using InvertedIndices
using PyCall

results = []

frames = [
    randn(Float32, 51, 51) .+ 100,
    randn(Float32, 101, 101) .+ 100,
    randn(Float32, 201, 201) .+ 100,
    randn(Float32, 301, 301) .+ 100,
    randn(Float32, 401, 401) .+ 100
]

# parameters
fwhm = 5


# julia benchmarks

for frame in frames
    N = length(frame)
    @info "Benchmarking - Julia" N=N

    time_elapsed = @belapsed detectionmap($frame, fwhm)
    
    @info "Benchmarking - Julia" N=N time=time_elapsed
    push!(results, (framework="ADI.jl", N=N, time=time_elapsed))
end

# python benchmarks

vip = pyimport("vip_hci")

for frame in frames[begin:3] # vip way too slow
    N = length(frame)
    @info "Benchmarking - Python" N=N

    time_elapsed = @belapsed vip.metrics.snrmap($frame, fwhm, verbose=false)
    
    @info "Benchmarking - Python" N=N time=time_elapsed
    push!(results, (framework="VIP", N=N, time=time_elapsed))
end

# load results file and update
path = joinpath(@__DIR__, "snrmap_benchmarks.csv")
df = isfile(path) ? CSV.File(path) |> DataFrame : DataFrame()
out = append!(DataFrame(results), df)
unique!(out, Not(:time)) |> CSV.write(path)
