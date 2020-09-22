using ProgressLogging
using Setfield

"""
    GreeDS(alg=TPCA(); threshold=0.0, progress=true)

Performs the greedy disk subtraction (GreeDS) algorithm.

This method is an iterative approach to standard ADI reduction which seeks to minimize over-subtraction by constraining the low-rank matrix approximation from `alg`.

For large data cubes the iteration can cause slowdowns, so a progress bar is provided using the [`ProgressLogging`](https://github.com/JunoLab/ProgressLogging.jl) API along with the `progress` keyword. It won't appear without a logging backend, such as [`TerminalLoggers`](https://github.com/c42f/TerminalLoggers.jl).

# Algorithms
Although the original paper explicitly uses PCA, we allow use of any linear ADI algorithm that is characterized by `ncomps`. By default, uses [`TPCA`](@ref).

# References
1. [Pairet et al. 2018](https://ui.adsabs.harvard.edu/abs/2018arXiv181201333P) "Reference-less algorithm for circumstellar disks imaging"
2. [Pairet et al. 2020](https://ui.adsabs.harvard.edu/abs/2020arXiv200805170P) "MAYONNAISE: a morphological components analysis pipeline for circumstellar disks and exoplanets imaging in the near infrared"
"""
struct GreeDS{ALG<:LinearAlgorithm} <: LinearAlgorithm
    alg::ALG
    threshold::Float64
    progress::Bool
end

GreeDS(alg=TPCA(); threshold=0.0, progress=true) = GreeDS(alg, threshold, progress)
GreeDS(ncomps::Int; kwargs...) = GreeDS(TPCA(ncomps); kwargs...)

function decompose(alg::GreeDS, cube, angles, cube_ref=cube; kwargs...)
    target = expand(cube)
    # get the number of components as a range from the underlying alg
    max_ncomps = isnothing(alg.alg.ncomps) ? size(target, 1) : alg.alg.ncomps
    # use the underlyhing algorithm with a lens for further processing
    _alg = alg.alg
    _alg = @set _alg.ncomps = 1
    reduced = _alg(target, angles, cube_ref)
    local design
    @progress for n in 1:max_ncomps
        resid = target .- expand_rotate(reduced, angles, alg.threshold; kwargs...)
        # use lens to update number of components
        _alg = @set _alg.ncomps = n
        # use the `resid` cube as the reference frames for the next reduction
        design = decompose(_alg, target, angles, resid; kwargs...)
        recon = reconstruct(_alg, design...)
        reduced = collapse!(target .- recon, angles; kwargs...)
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
