using ProgressLogging
using Setfield

"""
    Pairet(alg=PCA(); threshold=0.0)

Performs the bootstrapping algorithm defined in _Pairet et al. (2018)_.

This method is an iterative approach to standard ADI reduction which seeks to minimize over-subtraction by hand-crafting a reference library without any negative signal. This is accomplished with the internal `pairet_theta` function. The value for min-clipping is given by `threshold`.

# Algorithms
Although the original paper explicitly uses PCA, we allow use of any linear ADI algorithm that is characterized by `ncomps`. By default, uses [`PCA`](@ref).
* [`PCA`](@ref)
* `NMF` (not yet implemented)

# References
1. [Pairet et al. 2018 "Reference-less algorithm for circumstellar disks imaging"](https://ui.adsabs.harvard.edu/abs/2018arXiv181201333P)
"""
struct Pairet{ALG<:LinearAlgorithm} <: LinearAlgorithm
    alg::ALG
    threshold::Float64
end

Pairet(alg=PCA(); threshold=0.0) = Pairet(alg, threshold)

function decompose(alg::Pairet, cube, angles; kwargs...)
    # get the number of components as a range from the underlying alg
    max_ncomps = isnothing(alg.alg.ncomps) ? alg.alg.ncomps : size(cube, 1)
    # use the underlyhing algorithm with a lens for further processing
    _alg = alg.alg
    _alg = @set _alg.ncomps = 1
    reduced = _alg(cube, angles)
    local basis, weights
    @progress for n in 1:max_ncomps
        resid = cube .- pairet_theta(reduced, angles, alg.threshold; kwargs...)
        # use lens to update number of components
        _alg = @set _alg.ncomps = n
        # use the `resid` cube as the reference frames for the next reduction
        basis, weights = decompose(_alg, cube, angles, resid; kwargs...)
        recon = weights' * basis |> expand
        reduced = collapse!(cube .- recon, angles; kwargs...)
    end
    return basis, weights
end

"""
    pairet_theta(frame, angles, threshold; kwargs...)

Takes a frame, expands it into a cube, rotates it clockwise by `angles`, and min-clips at `threshold`. Keyword arguments will be passed to `HCIToolbox.derotate!`.
"""
function pairet_theta(frame, angles, threshold; kwargs...)
    N = length(angles)
    _frame = @. ifelse(frame > threshold, frame, threshold)
    cube = similar(frame, N, size(frame)...)
    for idx in axes(cube, 1)
        cube[idx, :, :] .= _frame
    end
    return derotate!(cube, -angles; kwargs...)
end
