using Statistics
using LinearAlgebra
using TSVD: tsvd


"""
    PCA(;ncomps=nothing, pratio=1) <: LinearAlgorithm

Use principal components analysis (PCA) to form a low-rank orthonormal basis of the input. Uses deterministic singular-value decomposition (SVD) to decompose data.

If `ncomps` is `nothing`, it will be set to the number of frames in the reference cube when processed.

# References
* [Soummer, Pueyo, and Larkin (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...755L..28S) "Detection and Characterization of Exoplanets and Disks Using Projections on Karhunen-Loève Eigenimages"

# Implements
* [`decompose`](@ref)
"""
@with_kw struct PCA <: LinearAlgorithm
    ncomps::Union{Int,Nothing} = nothing
    pratio::Float64 = 1.0
end

PCA(ncomps; kwargs...) = PCA(;ncomps=ncomps, kwargs...)

function decompose(alg::PCA, cube, angles, cube_ref=cube; kwargs...)
    @unpack ncomps, pratio = alg
    isnothing(ncomps) && (ncomps = size(cube, 1))
    ncomps > size(cube, 1) && error("ncomps ($ncomps) cannot be greater than the number of frames ($(size(cube, 1)))")

    # transform cube
    X = flatten(cube)
    X_ref = flatten(cube_ref)

    # fit SVD to get principal subspace of reference
    decomp = svd(X_ref)

    # get the minimum number comps to explain `pratio`
    pr = cumsum(decomp.S ./ sum(decomp.S))
    pr_n = findfirst(p -> p ≥ pratio, pr)
    nc = isnothing(pr_n) ? min(ncomps, size(decomp.Vt, 1)) : min(ncomps, pr_n, size(decomp.Vt, 1))

    nc < ncomps && @info "target pratio $pratio reached with only $nc components"
    # Get the principal components (principal subspace)
    P = decomp.Vt[1:nc, :]
    # reconstruct X using prinicipal subspace
    weights = P * X'

    return P, weights
end

"""
    TPCA(;ncomps=nothing) <: LinearAlgorithm

Perform principal components analysis (PCA) using truncated SVD (TSVD; provided by TSVD.jl) instead of deterministic SVD. This is often faster thant [`PCA`](@ref) but is non-determinstic. We find the differences unnoticable in practice.

If `ncomps` is `nothing`, it will be set to the number of frames in the reference cube when processed.

# Implements
* [`decompose`](@ref)

# See Also
* [`PCA`](@ref), [`TSVD.tsvd`](https://julialinearalgebra.github.io/TSVD.jl/latest/)
"""
@with_kw struct TPCA <: LinearAlgorithm
    ncomps::Union{Int,Nothing} = nothing
end

TPCA(ncomps) = TPCA(;ncomps=ncomps)

function decompose(alg::TPCA, cube, angles, cube_ref=cube; kwargs...)
    X = flatten(cube)
    X_ref = flatten(cube_ref)
    k = isnothing(alg.ncomps) ? size(cube, 1) : alg.ncomps
    _, _, P = tsvd(X_ref, k)
    A = P'
    w = A * X'
    return A, w
end
