using ProgressLogging
using Setfield

"""
    Pairet(alg=PCA(); threshold=0.0)

Performs the bootstrapping algorithm defined in _Pairet et al. (2018)_.

This method is an iterative approach to standard ADI reduction which seeks to minimize over-subtraction by hand-crafting a reference library without any negative signal. This is accomplished with the internal `pairet_theta` function. The value for min-clipping is given by `threshold`.

For large data cubes the iteration can cause slowdowns, so a progress bar is provided using the [`ProgressLogging`](https://github.com/JunoLab/ProgressLogging.jl) API. It won't appear without a logging backend, such as [`TerminalLoggers`](https://github.com/c42f/TerminalLoggers.jl).

# Algorithms
Although the original paper explicitly uses PCA, we allow use of any linear ADI algorithm that is characterized by `ncomps`. By default, uses [`PCA`](@ref).

# References
1. [Pairet et al. 2018](https://ui.adsabs.harvard.edu/abs/2018arXiv181201333P) "Reference-less algorithm for circumstellar disks imaging"
"""
struct Pairet{ALG<:LinearAlgorithm} <: LinearAlgorithm
    alg::ALG
    threshold::Float64
end

Pairet(alg=PCA(); threshold=0.0) = Pairet(alg, threshold)

function decompose(alg::Pairet, cube, angles; kwargs...)
    # get the number of components as a range from the underlying alg
    max_ncomps = isnothing(alg.alg.ncomps) ? size(cube, 1) : alg.alg.ncomps
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
        cube[idx, :, :] .= derotate(_frame, -angles[idx]; kwargs...)
    end
    return cube
end
