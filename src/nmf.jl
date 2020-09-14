using NMF: nnmf, solve!, CoordinateDescent, nndsvd

"""
    NMF(;ncomps=nothing) <: LinearAlgorithm
"""
struct NMF <: LinearAlgorithm
    ncomps::Int
end

function ADI.decompose(alg::NMF, cube, angles, cube_ref; kwargs...)
    X = flatten(cube)
    X_ref = flatten(cube_ref)
    W, H = nndsvd(X_ref, alg.ncomps)
    solve!(CoordinateDescent{eltype(X)}(maxiter=100, α=0), X, W, H)
    return H, W'
end

function ADI.decompose(alg::NMF, cube, angles; kwargs...)
    X = flatten(cube)
    res = nnmf(X, alg.ncomps; alg=:cd, init=:nndsvd, maxiter=100)
    return res.H, res.W'
end
