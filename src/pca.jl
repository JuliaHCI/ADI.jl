using Statistics
using LinearAlgebra

"""
    PCADesign(cube, [ref], angles; ncomps, pratio=1)

Create an object containing the decomposition of a cube (or reference) using principal component analysis (PCA) to form the approximate reconstruction of the systematics.

This will create a design matrix (the principal subspace) of the cube (or reference) truncated to either `ncomps` or until the prinicpal ratio is equal to `pratio` (whichever is fewer). As `ncomps` (or `pratio`) increase, more structure is removed from the cube, thus it is possible to over-subtract signal when choosing the size of the principal subpspace.

# See Alse
[`pca`](@ref), [`pairet`](@ref)
"""
struct PCADesign{T,C<:AbstractArray{T},M<:AbstractMatrix{T},V<:AbstractVector,N<:NamedTuple} <: ADIDesign{T, C, V}
    pcs::M
    weights::M
    S::C
    cube::C
    angles::V
    metadata::N
end

function Base.show(io::IO, d::PCADesign{T}) where T
    p, n = size(d.weights)
    print(io, "PCADesign{$T}(ncomps=$p, D=$n)")
    for (key, val) in pairs(d.metadata)
        print(io, "\n  $key: $val")
    end
    return nothing
end

design(d::PCADesign) = d.pcs'
weights(d::PCADesign) = d.weights
reconstruct(d::PCADesign) = reconstruct(d, d.cube)
reconstruct(d::PCADesign, X::AbstractMatrix) = (X * d.pcs') * d.pcs

PCADesign(cube::AbstractArray, angles::AbstractVector; kwargs...) = PCADesign(cube, cube, angles; kwargs...)

function PCADesign(cube::AbstractArray, ref::AbstractArray, angles::AbstractVector; ncomps, pratio = 1)
    # transform cube
    X_ref = flatten(ref)
    X = flatten(cube)

    # fit SVD to get principal subspace of reference
    decomp = svd(X_ref)

    # get the minimum number comps to explain `pratio`
    pr = cumsum(decomp.S ./ sum(decomp.S))
    pr_n = findfirst(p -> p ≥ pratio, pr)
    nc = pr_n === nothing ? min(ncomps, size(decomp.Vt, 1)) : min(ncomps, pr_n, size(decomp.Vt, 1))

    nc < ncomps && @info "target pratio $pratio reached with only $nc components"
    # Get the principal components (principal subspace)
    P = decomp.Vt[1:nc, :]
    # reconstruct X using prinicipal subspace
    weights = P * X'
    reconstructed = weights' * P |> expand
    S = cube .- reconstructed

    metadata = (;pratio=pr[nc])

    return PCADesign(promote(P, weights)..., promote(S, cube)..., normalize_par_angles(angles), metadata)
end

"""
    pca(cube, [ref], angles; ncomps, pratio=1, kwargs...)

Convenience function for [`PCADesign`](@ref) which returns the collapsed residual from the design. Any additional keyword arguments will be passed to [`reduce`](@ref).
"""
function pca(cube::AbstractArray, ref::AbstractArray, angles::AbstractVector; ncomps, pratio = 1, kwargs...)
    des = PCADesign(cube, ref, angles; ncomps=ncomps, pratio=pratio)
    result = reduce(des, kwargs...)
    return result
end

pca(cube::AbstractArray, angles::AbstractVector; kwargs...) =
    pca(cube, cube, angles; kwargs...)
