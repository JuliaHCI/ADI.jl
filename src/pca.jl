using TSVD: tsvd


"""
    PCA(ncomps; options...)
    PCA(;ncomps=nothing, options...)

Use principal components analysis (PCA) to form a low-rank orthonormal basis of the input. Uses deterministic singular-value decomposition (SVD) to decompose data.

If `ncomps` is `nothing`, the basis will not be truncated (i.e. `ncomps` is equal to the number of frames). `ncomps` can be set to `:noise` or `:pratio` to automatically choose the number of components using the residual frame noise or principal ratio, respectively. For more information, see the extended help.

# References
1. [Soummer, Pueyo, and Larkin (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...755L..28S) "Detection and Characterization of Exoplanets and Disks Using Projections on Karhunen-Loève Eigenimages"


# Extended help

## Optimizing `ncomps`

There are a few ways to optimize `ncomps` using the input data. Additional options for the optimization are listed below
1. `ncomps=:noise` - residual noise optimization
2. `ncomps=:pratio` - principal ratio optimization

### Residual noise optimization

This technique progressively increases `ncomps` at each step measuring the pixel-to-pixel noise (standard deviation) in the residual data. Iteration will stop when the noise is not improving beyond a threshold. This is suited for data with similar statistical characteristics, such as an annulus more so than a full-frame cube.
* `collapse=false` - if true, the temporal median of the residual data will be used for measuring the noise.
* `noise_error=1e-3` - the threshold for the minimal noise improvement looking back 2 iterations

### Principal ratio optimization

This technique chooses the number of components required to explain some ratio of the total variance in the data. This is known as the *principal ratio* or the *explained variance ratio*. The explained variance is measured by transforming the singular values of the SVD decomposition (`Λ = @. S^2 / (n - 1)`).
* `pratio=0.9` - the target principal ratio (between 0 and 1)
"""
@concrete struct PCA <: ADIAlgorithm
    ncomps
    opts
end
PCA(ncomps; options...) = PCA(ncomps, options)
PCA(; ncomps=nothing, options...) = PCA(ncomps, options)

function fit(alg::PCA, data::AbstractMatrix; ref=data, kwargs...)
    # get number of components (using dispatch for symbolic args)
    k = get_ncomps(alg.ncomps, ref; alg.opts...)
    # fit SVD to get principal subspace of reference
    decomp = svd!(collect(ref))
    # Get the principal components (principal subspace) and weights
    P = decomp.Vt[1:k, :]
    weights = data * P'
    return LinearDesign(P, weights)
end

# get ncomps using given value or num frames, whichever is smaller
get_ncomps(n::Int, data; kwargs...) = min(n, size(data, 1))
get_ncomps(::Nothing, data; kwargs...) = size(data, 1)

# get ncomps using automatic methods
function get_ncomps(s::Symbol, data; kwargs...)
    if s === :noise
        noise_decay_ncomps(data; kwargs...)
    elseif s === :pratio
        pratio_ncomps(data; kwargs...)
    else
        error("Invalid `ncomps`. Did you mean :noise or :pratio?")
    end
end

function noise_decay_ncomps(data; collapse=false, noise_error=1e-3)
    if collapse
        μ = mean(data; dims=1)
        σ2 = var(data; dims=1, mean=μ)
        X = @. (data - μ) / σ2
    else
        X = data .- mean(data, dims=1)
    end
    P = svd(X).Vt
    tmpr = similar(data)
    τ1 = τ2 = 0
    @progress name="Optimizing ncomps using residual noise" for ncomp in axes(data, 1)
        Pv = @view P[1:ncomp, :]
        tmpr .= X * (I - Pv'Pv)
        # calculate noise (standard deviation) optionally collapsing
        noise = collapse ? std(median(tmpr, dims=1)) : std(tmpr)
        # test if we've reached the noise decay tolerance
        if ncomp > firstindex(data, 1) + 2
            px_noise_decay = τ2 - noise
            @debug noise_decay=px_noise_decay noise=noise
            if px_noise_decay < noise_error
                @info "noise threshold reached with $ncomp components"
                return ncomp
            end
        end
        # update recursion variables
        τ2, τ1 = τ1, noise
    end
    return lastindex(data, 1)
end

function pratio_ncomps(data; pratio=0.9)
    @debug "Choosing ncomps required to explain $(pratio*100)% of data's temporal variance"
    X = data .- mean(data, dims=1)
    Λ = svd!(X).S
    n = length(Λ)
    exp_var = @. Λ^2 / (n - 1)
    ratio_cumsum = cumsum(exp_var ./ sum(exp_var))
    k = last(searchsorted(ratio_cumsum, pratio))
    @info "$(pratio*100)% of variance explained with $k components"
    return k
end

"""
    TPCA(;ncomps=nothing, options...)
    TPCA(ncomps; options...)

Perform principal components analysis (PCA) using truncated SVD (TSVD; provided by TSVD.jl) instead of deterministic SVD. This is often faster than [`PCA`](@ref) but is non-deterministic, so the results may be different. See the [`PCA`](@ref) docstring for more information about the options.

# See Also
* [`PCA`](@ref), [`TSVD.tsvd`](https://github.com/JuliaLinearAlgebra/TSVD.jl)
"""
@concrete struct TPCA <: ADIAlgorithm
    ncomps
    opts
end
TPCA(ncomps; options...) = TPCA(ncomps, options)
TPCA(; ncomps=nothing, options...) = TPCA(ncomps, options)

function fit(alg::TPCA, data::AbstractMatrix; ref=data, kwargs...)
    # get number of components (using dispatch for symbolic args)
    # TODO using automatic methods is not valid for TPCA using different SVD methods
    k = get_ncomps(alg.ncomps, ref; alg.opts...)
    # fit SVD to get principal subspace of reference
    U, Σ, V = tsvd(collect(ref), k)
    # Get the principal components (principal subspace) and weights
    P = V'
    weights = data * P'
    return LinearDesign(P, weights)
end
