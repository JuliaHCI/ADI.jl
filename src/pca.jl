using Statistics
using LinearAlgebra


"""
    PCA(;ncomps=nothing, pratio=1)

If `ncomps` is `nothing`, it will be set to the number of frames in the reference cube when processed.

# References
* [Soummer, Pueyo, and Larkin (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...755L..28S) "Detection and Characterization of Exoplanets and Disks Using Projections on Karhunen-Loève Eigenimages"
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
    X = ndims(cube) === 2 ? cube : flatten(cube)
    X_ref = ndims(cube_ref) === 2 ? cube_ref : flatten(cube_ref)

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
