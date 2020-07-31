using ProgressLogging

"""
    pairet([alg=pca], cube, angles; ncomps, threshold=0, kwargs...)

Performs the bootstrapping algorithm defined in Pairet et al. 2018.

This method is an iterative approach to standard ADI reduction which seeks to minimze over-subtracting signal in the low-rank approximation of `cube`.

`ncomps` can be an integer, which will iterate over `1:ncomps`, otherwise it can be any sub-type of `AbstractRange{<:Int}`. As part of the Pairet algorithm, the low-rank approximation at each iteration will be min-clipped at `threshold`.        

Any extra keyword arguments will be passed to `reduce(::ADIDesign)`.

# Algorithms
Although the original paper explicitly uses PCA, we allow use of any ADI algorithm that is characterized by `ncomps`. By default, uses [`pca`](@ref).
* [`pca`](@ref)
* `nmf` (not yet implemented)

### References
1. [Pairet et al. 2018 "Reference-less algorithm for circumstellar disks imaging"](https://ui.adsabs.harvard.edu/abs/2018arXiv181201333P)
"""
function pairet(alg, cube::AbstractArray{T}, angles::AbstractVector; ncomps, threshold=zero(T), kwargs...) where T
    design = alg(cube, angles; ncomps=1)
    red = reduce(design; kwargs...)

    @progress for n in _get_range(ncomps)
        resid = cube .- _pairet_theta(red, design.angles, threshold; kwargs...)
        design = alg(cube, resid, angles; ncomps=n)
        red .= reduce(design; kwargs...)
    end
    return design
end

pairet(cube::AbstractArray, angles::AbstractVector; kwargs...) = pairet(PCADesign, cube, angles; kwargs...)

_get_range(n::Integer) = 1:n
_get_range(n::AbstractRange{<:Integer}) = n

# takes a frame, expands it into a cube, rotates it clockwise by angles, 
# and min-clips at threshold
function _pairet_theta(frame, angles, threshold; kwargs...)
    N = length(angles)
    _frame = @. ifelse(frame > threshold, frame, threshold)
    cube = similar(frame, N, size(frame)...)
    for idx in axes(cube, 1)
        cube[idx, :, :] .= _frame
    end
    return derotate!(cube, -angles; kwargs...)
end
