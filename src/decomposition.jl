import MultivariateStats
using Statistics
using LinearAlgebra: I
using NMF: nnmf

const mvs = MultivariateStats




# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------

"""
    NMF

Use non-negative matrix factorization (NMF) to reduce data cube.

Uses [`NMF.nnmf`](https://github.com/JuliaStats/NMF.jl) for decomposition.

# Arguments
* `ncomps::Int` - The number of components to keep. Cannot be larger than the number of frames in the input cube (default).

# Examples

```jldoctest
julia> cube = ones(30, 100, 100);

julia> design(NMF, cube)
(A = [0.18699723925131992 0.0017227311439991433 … 0.25104857267542474 0.007639437531209008; 0.17892712604685682 0.0017227311439991433 … 0.5290547113213493 0.0039969810485154245; … ; 0.1791108604483525 0.3333726515006622 … 0.009151289339269827 0.0039969810485154245; 0.18321575873721374 0.3333726515006622 … 0.009151289339269827 2.90943813485493e-13], w = [4.706954006158167 4.706954006158167 … 4.706954006158167 4.706954006158167; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0], S = [0.9999999998346566 0.9999999998346566 … 0.9999999998346566 0.9999999998346566; 1.0000000001363376 1.0000000001363376 … 1.0000000001363376 1.0000000001363376; … ; 1.0000000001296059 1.0000000001296059 … 1.0000000001296059 1.0000000001296059; 0.9999999999760179 0.9999999999760179 … 0.9999999999760179 0.9999999999760179])

julia> design(NMF, cube, 5)
(A = [0.18257421287124714 0.005547894341040443 … 0.09696771072251913 0.15665986349016803; 0.1825741813281496 0.0055478943410398615 … 0.012760078325461317 0.25756981304356596; … ; 0.18257415421552403 0.33337265150066164 … 0.006436229867279868 0.32699785187779695; 0.18257415421552403 0.33337265150066164 … 0.006436229867279868 0.32699785187779695], w = [5.47722321153941 5.47722321153941 … 5.47722321153941 5.47722321153941; 1.0094317927007388e-7 1.0094317927007388e-7 … 1.0094317927007388e-7 1.0094317927007388e-7; … ; 6.241103130392369e-8 6.241103130392369e-8 … 6.241103130392369e-8 6.241103130392369e-8; 1.7433266010891995e-6 1.7433266010891995e-6 … 1.7433266010891995e-6 1.7433266010891995e-6], S = [0.9999999999004122 0.9999999999004122 … 0.9999999999004122 0.9999999999004122; 1.0000000000166112 1.0000000000166112 … 1.0000000000166112 1.0000000000166112; … ; 1.0000000001164895 1.0000000001164895 … 1.0000000001164895 1.0000000001164895; 1.0000000001164895 1.0000000001164895 … 1.0000000001164895 1.0000000001164895])

```
"""
struct NMF <: HCIAlgorithm end

function design(::Type{<:NMF}, cube::AbstractArray{T,3}, ncomps::Integer = size(cube, 1); kwargs...) where T
    flat_cube = flatten(cube)

    nmf = nnmf(flat_cube, ncomps; kwargs...)
    nmf.converged || @warn "NMF did not converge, try changing `alg`, `maxiter` or `tol` as keyword wargs."

    A = nmf.W
    weights = nmf.H
    reconstructed = nmf.W * nmf.H
    return (A = A, w = weights, S = reconstructed)
end

# ------------------------------------------------------------------------------

"""
    Pairet{<:Union{PCA, NMF}}
"""
struct Pairet{T <: HCIAlgorithm} <: HCIAlgorithm end


function design(::Type{Pairet{<:D}}, cube::AbstractArray{T,3}, ncomps::Integer = size(cube, 1); pca_kwargs...) where {D <: Union{PCA,NMF},T}
    S = []
    reduced = []

    X = flatten(cube)

initial_pca = mvs.fit(mvs.PCA, X; maxoutdim = 1)
    # TODO
end

# ------------------------------------------------------------------------------

"""
    Median

Design using the median of the cube

# Examples
```jldoctest
julia> cube = ones(30, 100, 100);

julia> design(Median, cube)
(A = [1.0 1.0 … 1.0 1.0], w = LinearAlgebra.UniformScaling{Bool}(true), S = [1.0 1.0 … 1.0 1.0])
```

# See Also
[`Mean`](@ref)
"""
struct Median <: HCIAlgorithm end

function design(::Type{<:Median}, cube::AbstractArray{T,3}) where T
    out = median(flatten(cube), dims = 1)
    weights = I
    return (A = out, w = weights, S = out)
end


"""
    Mean

Design using the mean of the cube

# Examples
```jldoctest
julia> cube = ones(30, 100, 100);

julia> design(Mean, cube)
(A = [1.0 1.0 … 1.0 1.0], w = LinearAlgebra.UniformScaling{Bool}(true), S = [1.0 1.0 … 1.0 1.0])
```

# See Also
[`Median`](@ref)
"""
struct Mean <: HCIAlgorithm end

function design(::Type{<:Mean}, cube::AbstractArray{T,3}) where T
    out = mean(flatten(cube), dims = 1)
    weights = I
    return (A = out, w = weights, S = out)
end


# ------------------------------------------------------------------------------

