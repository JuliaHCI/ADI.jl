using FillArrays
using Setfield

const PCALIKE = Union{PCA,NMF}

"""
    GreeDS(alg=PCA(); threshold=0)
    GreeDS(ncomps; threshold=0, options...)

Performs the greedy disk subtraction (GreeDS) algorithm.

This method is an iterative approach to standard ADI reduction which seeks to minimize over-subtraction by constraining the low-rank matrix approximation from `alg`. By default, uses [`PCA`](@ref). If `ncomps` or other PCA options are provided, they will be passed to the constructor.

!!! note
    The GreeDS algorithm requires fully reconstructing a cube at each iteration, which requires knowing the geometry of the input (full-frame, annulus, etc.) and the corresponding parallactic angles. These angles must be passed as a keyword argument `angles`. In the case of reducing data, e.g. `GreeDS()(cube, angles)` the angles will be passed automatically. It is important to clarify, *these angles should correspond to the reference data in the case of RDI*, e.g. `GreeDS()(cube, angles; ref=ref_cube, angles=ref_angles)`

# Algorithms
The following algorithms work natively with GreeDS: [`PCA`](@ref) and [`NMF`](@ref)

# References
1. [Pairet et al. 2018](https://ui.adsabs.harvard.edu/abs/2018arXiv181201333P) "Reference-less algorithm for circumstellar disks imaging"
2. [Pairet et al. 2020](https://ui.adsabs.harvard.edu/abs/2020arXiv200805170P) "MAYONNAISE: a morphological components analysis pipeline for circumstellar disks and exoplanets imaging in the near infrared"
"""
@concrete struct GreeDS{AT<:PCALIKE} <: ADIAlgorithm
    kernel::AT
    threshold
end
GreeDS(alg=PCA(); threshold=0) = GreeDS(alg, threshold)
GreeDS(ncomps::Int; threshold=0, kwargs...) = GreeDS(PCA(ncomps; kwargs...), threshold=threshold)

function fit(alg::GreeDS, data::AbstractArray{T,3}; angles, ref::AbstractArray{S,3}=data, kwargs...) where {T,S}
    target = flatten(ref)
    # get the number of components as a range from the underlying alg
    max_ncomps = get_ncomps(alg.kernel.ncomps, target)
    # use the underlyhing algorithm with a lens for further processing
    f = alg.kernel
    f = @set f.ncomps = 1
    design = fit(f, target)
    R = expand(target .- reconstruct(design))
    reduced = collapse!(R, angles)
    @progress name="GreeDS" for n in 1:max_ncomps
        resid = ref .- expand_rotate(reduced, angles, alg.threshold)
        # use lens to update number of components
        f = @set f.ncomps = n
        # use the `resid` cube as the reference frames for the next reduction
        design = fit(f, target; ref=flatten(resid))
        R = expand(target .- reconstruct(design))
        reduced = collapse!(R, angles)
    end
    # RDI not defined in Pairet 18,20; project onto components
    if ref !== data
        A = design.basis
        weights = transpose(A) * flatten(data)
        return LinearDesign(A, weights)
    end
    return design
end

function fit(alg::GreeDS, data::AnnulusView; angles, ref::AnnulusView=data, kwargs...)
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
    @progress name="GreeDS" for n in 1:max_ncomps
        tmpAnn .= ref .- expand_rotate(reduced, angles, alg.threshold)
        resid = tmpAnn()
        # use lens to update number of components
        f = @set f.ncomps = n
        # use the `resid` cube as the reference frames for the next reduction
        design = fit(f, target; ref=resid)
        R = inverse(data, target .- reconstruct(design))
        reduced = collapse!(R, angles)
    end
    if ref !== data
        A = design.basis
        weights = transpose(A) * data()
        return LinearDesign(A, weights)
    end
    return design
end

"""
    expand_rotate(frame, angles, threshold; kwargs...)

Takes a frame, expands it into a cube, rotates it clockwise by `angles`, and min-clips at `threshold`. Keyword arguments will be passed to `HCIToolbox.derotate`.
"""
function expand_rotate(frame, angles, threshold; kwargs...)
    N = length(angles)
    _frame = max.(frame, threshold)
    cube = similar(frame, size(frame)..., N)
    Threads.@threads for idx in axes(cube, 3)
        cube[:, :, idx] .= derotate(_frame, -angles[idx]; kwargs...)
    end
    return cube
end
