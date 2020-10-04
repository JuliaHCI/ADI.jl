using NMF: nnmf, solve!, CoordinateDescent, nndsvd

"""
    NMF(;ncomps=nothing) <: LinearAlgorithm
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
