module ADIX

using LinearAlgebra
using NMF: nndsvd, solve!, CoordinateDescent
using Parameters
using ProgressLogging
using Reexport
using Setfield
using Statistics
using StructArrays

@reexport using HCIToolbox

export reconstruct,
       process,
       Classic,
       PCA,
       NMF,
       GreeDS


#####################
###   INTERFACE   ###
#####################


abstract type ADIAlgorithm end

"""
    PCA(ncomps=nothing; options...)
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
struct PCA <: ADIAlgorithm
    ncomps
    opts
end
PCA(ncomps; options...) = PCA(ncomps, options)
PCA(; ncomps=nothing, options...) = PCA(ncomps, options)

struct NMF <: ADIAlgorithm
    ncomps
end
NMF(; ncomps=nothing) = NMF(ncomps)

struct Classic <: ADIAlgorithm
    method
end
Classic(;method=median) = Classic(method)


"""
    GreeDS(alg=PCA(); threshold=0.0)
    GreeDS(ncomps; threshold=0.0, options...)

Performs the greedy disk subtraction (GreeDS) algorithm.

This method is an iterative approach to standard ADI reduction which seeks to minimize over-subtraction by constraining the low-rank matrix approximation from `alg`. By default, uses [`PCA`](@ref). If `ncomps` or other PCA options are provided, they will be passed to the constructor.

For large data cubes the iteration can cause slowdowns, so a progress bar is provided using the [`ProgressLogging`](https://github.com/JunoLab/ProgressLogging.jl) API along with the `progress` keyword. It won't appear without a logging backend, such as [`TerminalLoggers`](https://github.com/c42f/TerminalLoggers.jl).

# Algorithms
Originally multiple algorithms were supported, but currently only [`PCA`](@ref) works properly with the GreeDS algorithm. 

# References
1. [Pairet et al. 2018](https://ui.adsabs.harvard.edu/abs/2018arXiv181201333P) "Reference-less algorithm for circumstellar disks imaging"
2. [Pairet et al. 2020](https://ui.adsabs.harvard.edu/abs/2020arXiv200805170P) "MAYONNAISE: a morphological components analysis pipeline for circumstellar disks and exoplanets imaging in the near infrared"
"""
struct GreeDS{AT} <: ADIAlgorithm
    kernel::AT
    threshold
end
GreeDS(alg=PCA(); threshold=0.0) = GreeDS(alg, threshold)
GreeDS(ncomps::Int; threshold=0.0, kwargs...) = GreeDS(PCA(ncomps; kwargs...), threshold=threshold)

function fit(alg::ADIAlgorithm, cube::AbstractArray{T,3}; kwargs...) where T
    data = flatten(cube)
    if :ref in keys(kwargs)
        kwargs[:ref] isa AbstractArray{<:Any,3} || error("reference data geometry does not match target data")
        ref_data = flatten(kwargs[:ref])
        return fit(alg, data; kwargs..., ref=ref_data)
    end
    return fit(alg, data; kwargs...)
end
function fit(alg::ADIAlgorithm, cube::AnnulusView; kwargs...)
    data = cube()
    if :ref in keys(kwargs)
        kwargs[:ref] isa AnnulusView || error("reference data geometry does not match target data")
        ref_data = kwargs[:ref]()
        return fit(alg, data; kwargs..., ref=ref_data)
    end
    return fit(alg, data; kwargs...)
end

function fit(alg::ADIAlgorithm, cube::MultiAnnulusView; kwargs...)
    if :ref in keys(kwargs)
        kwargs[:ref] isa MultiAnnulusView || error("reference data geometry does not match target data")
        anns = eachannulus(cube)
        ref_anns = eachannulus(kwargs[:ref])
        return StructArray(fit(alg, ann; kwargs..., ref=ref_ann) for (ann, ref_ann) in zip(anns, ref_anns))
    else
        return StructArray(fit(alg, ann; kwargs...) for ann in eachannulus(cube))
    end
end

function reconstruct(alg::ADIAlgorithm, cube; kwargs...)
    design = fit(alg, cube; kwargs...)
    S = reconstruct(design)
    return expand_geometry(cube, S)
end

expand_geometry(::AbstractArray{T,3}, arr) where {T} = expand(arr)
expand_geometry(cube::AnnulusView, arr) = inverse(cube, arr)
expand_geometry(cube::MultiAnnulusView, arrs) = inverse(cube, arrs)

function process(alg::ADIAlgorithm, cube; kwargs...)
    S = reconstruct(alg, cube; kwargs...)
    return cube .- S
end

function process(alg::NMF, cube; kwargs...)
    if any(<(0), cube)
        target = copy(cube)
        target .-= minimum(target)
    else
        target = cube
    end
    S = reconstruct(alg, target; kwargs...)
    return target .- S
end

function (alg::ADIAlgorithm)(cube, angles; method=:deweight, kwargs...)
    return collapse!(process(alg, cube; kwargs...), angles, method=method)
end

function (alg::GreeDS)(cube, angles; method=:deweight, kwargs...)
    return collapse!(process(alg, cube; angles=angles, kwargs...), angles, method=method)
end

###########################
###   DATA STRUCTURES   ###
###########################

abstract type AbstractDesign end
abstract type LinearDesign <: AbstractDesign end

struct PCADesign{AT,WT} <: LinearDesign
    ncomps::Int
    components::AT
    weights::WT
end

struct NMFDesign{AT,WT} <: LinearDesign
    components::AT
    weights::WT
end

struct ClassicDesign{FT} <: AbstractDesign
    n::Int # number of frames in original data
    frame::FT
end

design(des::ClassicDesign) = des.frame
design(des::PCADesign) = (des.components, des.weights)
design(des::NMFDesign) = (des.components, des.weights)

function reconstruct(des::LinearDesign)
    A, weights = design(des)
    return weights * A
end
reconstruct(des::ClassicDesign) = repeat(des.frame, des.n, 1)
reconstruct(designs::AbstractVector{<:AbstractDesign}) = map(reconstruct, designs)

######################
###   ALGORITHMS   ###
######################

function fit(alg::PCA, data::AbstractMatrix; ref=data, kwargs...)
    # get number of components (using dispatch for symbolic args)
    k = get_ncomps(alg.ncomps, ref; alg.opts...)
    # fit SVD to get principal subspace of reference
    decomp = svd(ref)
    # Get the principal components (principal subspace) and weights
    P = decomp.Vt[begin:k, :]
    weights = data * P'
    return PCADesign(k, P, weights)
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
    @progress "Optimizing ncomps using residual noise" for ncomp in axes(data, 1)
        Pv = @view P[begin:ncomp, :]
        tmpr .= X * (I - Pv'Pv)
        # calculate noise (standard deviation) optionally collapsing
        noise = collapse ? std(median(tmpr, dims=1)) : std(tmpr)
        # test if we've reached the noise decay tolerance
        if ncomp > firstindex(data, 1) + 2
            px_noise_decay = τ2 - noise
            @debug noise_decay=px_noise_decay noise=noise
            px_noise_decay < noise_error && return ncomp
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
    return last(searchsorted(ratio_cumsum, pratio))
end

function fit(alg::Classic, data::AbstractMatrix; ref=data)
    n = size(data, 1)
    return ClassicDesign(n, alg.method(ref, dims=1))
end

function fit(alg::NMF, data::AbstractMatrix{T}) where T
    if any(<(0), data)
        @warn "Negative values encountered in input. Make sure to rescale your data"
    end
    k = min(alg.ncomps, size(data, 1))

    W, H = nndsvd(data, k)
    solve!(CoordinateDescent{T}(), data, W, H)
    return NMFDesign(H, W)
end



function fit(alg::GreeDS{<:PCA}, data::AbstractArray{T,3}; angles, ref=data) where T
    target = flatten(ref)
    # get the number of components as a range from the underlying alg
    max_ncomps = get_ncomps(alg.kernel.ncomps, target)
    # use the underlyhing algorithm with a lens for further processing
    f = alg.kernel
    f = @set f.ncomps = 1
    design = fit(f, target)
    R = expand(target .- reconstruct(design))
    reduced = collapse!(R, angles)
    @progress "GreeDS" for n in 1:max_ncomps
        resid = data .- expand_rotate(reduced, angles, alg.threshold)
        # use lens to update number of components
        f = @set f.ncomps = n
        # use the `resid` cube as the reference frames for the next reduction
        design = fit(f, target; ref=flatten(resid))
        R = expand(target .- reconstruct(design))
        reduced = collapse!(R, angles)
    end
    if ref !== data
        A = design.components
        design.weights .= flatten(data) * A'A
    end
    return design
end

function fit(alg::GreeDS{<:PCA}, data::AnnulusView; angles, ref=data)
    target = data()
    tmpAnn = copy(data)
    # get the number of components as a range from the underlying alg
    max_ncomps = get_ncomps(alg.kernel.ncomps, target)
    # use the underlyhing algorithm with a lens for further processing
    f = alg.kernel
    f = @set f.ncomps = 1
    design = fit(f, target)
    R = inverse(data, target .- reconstruct(design))
    reduced = collapse!(R, angles)
    @progress "GreeDS" for n in 1:max_ncomps
        tmpAnn.parent .= data .- expand_rotate(reduced, angles, alg.threshold)
        resid = tmpAnn()
        # use lens to update number of components
        f = @set f.ncomps = n
        # use the `resid` cube as the reference frames for the next reduction
        design = fit(f, target; ref=resid)
        design.weights .= target * design.components'
        R = inverse(data, target .- reconstruct(design))
        reduced = collapse!(R, angles)
    end
    return design
end


"""
    expand_rotate(frame, angles, threshold; kwargs...)

Takes a frame, expands it into a cube, rotates it clockwise by `angles`, and min-clips at `threshold`. Keyword arguments will be passed to `HCIToolbox.derotate!`.
"""
function expand_rotate(frame, angles, threshold; kwargs...)
    N = length(angles)
    _frame = @. ifelse(frame > threshold, frame, threshold)
    cube = similar(frame, N, size(frame)...)
    Threads.@threads for idx in axes(cube, 1)
        cube[idx, :, :] .= derotate(_frame, -angles[idx]; kwargs...)
    end
    return cube
end


end # module
