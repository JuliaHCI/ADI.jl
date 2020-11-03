using NMF: nnmf, solve!, CoordinateDescent, nndsvd

"""
    NMF(;ncomps=nothing) <: LinearAlgorithm
    NMF(ncomps)

Use non-negative matrix factorization (NMF, NNMF) to form a non-negative low-rank orthonormal basis of the input. The implementation of the underlying NMF is provided by [NMF.jl](https://github.com/JuliaStats/NMF.jl). The implementation uses a non-negative SVD for initialization and a coordinate-descent solver to fit.

If `ncomps` is `nothing`, it will be set to the number of frames in the reference cube when processed.

# References
* [Ren et al. 2018](http://adsabs.harvard.edu/abs/2018ApJ...852..104R) Non-negative Matrix Factorization: Robust Extraction of Extended Structures
"""
@with_kw struct NMF <: LinearAlgorithm
    ncomps::Union{Int,Nothing} = nothing
end

function ADI.decompose(alg::NMF, cube, angles, cube_ref; kwargs...)
    X = flatten(cube)
    X_ref = flatten(cube_ref)

    k = isnothing(alg.ncomps) ? size(cube, 1) : alg.ncomps
    k > size(cube, 1) && error("ncomps ($k) cannot be greater than the number of frames ($(size(cube, 1)))")
    # manually initialize W and H
    _, H = nndsvd(X_ref, k)
    W = X * H'
    solve!(CoordinateDescent{eltype(X)}(maxiter=100, α=0), X, W, H)
    return H, W
end

function ADI.decompose(alg::NMF, cube, angles; kwargs...)
    X = flatten(cube)
    
    k = isnothing(alg.ncomps) ? size(cube, 1) : alg.ncomps
    k > size(cube, 1) && error("ncomps ($k) cannot be greater than the number of frames ($(size(cube, 1)))")

    res = nnmf(X, k; alg=:cd, init=:nndsvd, maxiter=100)
    return res.H, res.W
end
