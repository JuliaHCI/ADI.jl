using Statistics
using LinearAlgebra

struct PCADesign{T<:AbstractArray,M<:AbstractMatrix,V<:AbstractVector} <: ADIDesign{T, V}
    A::M
    w::M
    reconstructed::T
    S::T
    angles::V
end

function Base.show(io::IO, d::PCADesign{T}) where T
    p, n = size(d.w)
    print("PCADesign{$T}(ncomps=$p, D=$n)")
end

reconstruct(d::PCADesign, X::AbstractMatrix) = (X * d.A') * d.A

"""
    pca(cube, [ref], angles, ncomps; pratio=0.99)

Use principal component analysis (PCA) to reduce data cube.
"""
pca(cube::AbstractArray, angles::AbstractVector; ncomps, pratio = 1) = pca(cube, cube, angles; ncomps=ncomps, pratio=pratio)

function pca(cube::AbstractArray, ref::AbstractArray, angles::AbstractVector; ncomps, pratio = 1)
    # transform cube
    X_ref = flatten(ref)
    X = flatten(cube)

    # fit SVD to get principal subspace
    U, = svd(X_ref')
    P = U[:, 1:ncomps]'

    # reconstruct X using prinicipal subspace
    weights = P * X'
    reconstructed = weights' * P |> expand

    S = cube .- reconstructed
    return PCADesign(collect(P), weights, reconstructed, S, normalize_par_angles(angles))
end
