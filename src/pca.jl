using Statistics
using LinearAlgebra

"""
    PCADesign(A, w, reconstructed, S, angles, pratio)

A container for PCA-based ADI algorithm output.

`A` is the principal subspace, `w` are the weights used to reconstruct to the target cube. The residual from reconstruction is stored in `S` as a cube. The parallactic-angles are stored in `angles`. Finally, the principal ratio (the ratio of explained variance per component to the total explained variance) is stored in `pratio`. 

# See Alse
[`pca`](@ref), [`pairet`](@ref)
"""
struct PCADesign{T<:AbstractArray,M<:AbstractMatrix,V<:AbstractVector,F<:AbstractFloat} <: ADIDesign{T, V}
    A::M
    w::M
    S::T
    angles::V
    pratio::F
end

function Base.show(io::IO, d::PCADesign{T}) where T
    p, n = size(d.w)
    print("PCADesign{$T}(ncomps=$p, D=$n, pratio=$(d.pratio))")
    return nothing
end

reconstruct(d::PCADesign, X::AbstractMatrix) = (X * d.A') * d.A

"""
    pca(cube, [ref], angles; ncomps, pratio=1)::PCADesign

Decomposes a cube (or reference) using principal component analysis (PCA) to form the approximate reconstruction of the systematic noise.

This will create a design matrix (the principal subspace) of the cube (or reference) truncated to either `ncomps` or until the prinicpal ratio is equal to `pratio` (whichever is fewer). As `ncomps` (or `pratio`) increase, more structure is removed from the cube, thus it is possible to over-subtract signal when choosing the size of the principal subpspace.
"""
pca(cube::AbstractArray, angles::AbstractVector; ncomps, pratio = 1) = pca(cube, cube, angles; ncomps=ncomps, pratio=pratio)

function pca(cube::AbstractArray, ref::AbstractArray, angles::AbstractVector; ncomps, pratio = 1)
    # transform cube
    X_ref = flatten(ref)
    X = flatten(cube)

    # fit SVD to get principal subspace of reference
    U, S, = svd(X_ref')

    # get the minimum number comps to explain `pratio`
    pr = cumsum(S ./ sum(S))
    pr_n = findfirst(p -> p ≥ pratio, pr)
    nc = pr_n === nothing ? min(ncomps, size(U, 2)) : min(ncomps, pr_n, size(U, 2))

    # Get the principal components (principal subspace)
    P = U[:, 1:nc]' |> collect
    # reconstruct X using prinicipal subspace
    weights = P * X'
    reconstructed = weights' * P |> expand
    S = cube .- reconstructed

    return PCADesign(promote(P, weights)..., S, normalize_par_angles(angles), pr[nc])
end
