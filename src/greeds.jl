using ProgressLogging
using Setfield

"""
    GreeDS(alg=PCA(); threshold=0.0)
    GreeDS(ncomps; threshold=0.0, options...)

Performs the greedy disk subtraction (GreeDS) algorithm.

This method is an iterative approach to standard ADI reduction which seeks to minimize over-subtraction by constraining the low-rank matrix approximation from `alg`. By default, uses [`PCA`](@ref). If `ncomps` or other PCA options are provided, they will be passed to the constructor.

For large data cubes the iteration can cause slowdowns, so a progress bar is provided using the [`ProgressLogging`](https://github.com/JunoLab/ProgressLogging.jl) API along with the `progress` keyword. It won't appear without a logging backend, such as [`TerminalLoggers`](https://github.com/c42f/TerminalLoggers.jl).

!!! note
    The GreeDS algorithm requires fully reconstructing a cube at each iteration, which requires knowing the geometry of the input (full-frame, annulus, etc.) and the corresponding parallactic angles. These angles must be passed as a keyword argument `angles`. In the case of reducing data, e.g. `GreeDS()(cube, angles)` the angles will be passed automatically.

# Algorithms
Originally multiple algorithms were supported, but currently only [`PCA`](@ref) and [`TPCA`](@ref) work properly with the GreeDS algorithm.

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

# make sure angles get passed automatically if doing full processing
function (alg::GreeDS)(cube, angles; method=:deweight, kwargs...)
    return collapse!(subtract(alg, cube; angles=angles, kwargs...), angles, method=method)
end

function fit(alg::GreeDS{<:Union{PCA,TPCA}}, data::AbstractArray{T,3}; angles, ref=data) where T
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
        weights .= flatten(data) * A'A
        return PCADesign(max_ncomps, A, weights)
    end
    return design
end

function fit(alg::GreeDS{<:Union{PCA,TPCA}}, data::AnnulusView; angles, ref=data)
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
        R = inverse(data, target .- reconstruct(design))
        reduced = collapse!(R, angles)
    end
    if ref !== data
        A = design.components
        weights = data() * A'A
        return PCADesign(max_ncomps, A, weights)
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
