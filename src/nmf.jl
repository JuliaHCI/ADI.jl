using NMF: nnmf, solve!, CoordinateDescent, nndsvd

"""
    NMF(;ncomps=nothing) <: LinearAlgorithm
    NMF(ncomps)

Use non-negative matrix factorization (NMF, NNMF) to form a non-negative, low-rank, and orthonormal basis of the input. The implementation of the underlying NMF is provided by [NMF.jl](https://github.com/JuliaStats/NMF.jl). The implementation uses a non-negative SVD for initialization and a coordinate-descent solver to fit.

If `ncomps` is `nothing`, it will be set to the number of frames in the reference cube when processed.

!!! warning "Non-negativity constraint"
    NMF is not designed to fit negative values. This algorithm will warn you (but will not error) if a target or reference cube contains negative values. The output may seem reasonable, but it is not well-defined with respect to the NMF algorithm. To overcome this, rescaling the data by its minimum before processing is necessary
    ```julia
    target = cube .- minimum(cube)
    S = reconstruct(NMF(), target, angles)
    ```
    When doing full-frame reduction (e.g. `NMF()(cube, angles)`) *this is handled automatically*, so this constraint only applies to the lower-level API and methods which rely on those, such as [`GreeDS`](@ref). **In general, if you see warnings, heed them.**

# Implements
* [`decompose`](@ref)

# References
1. [Ren et al. 2018](http://adsabs.harvard.edu/abs/2018ApJ...852..104R) Non-negative Matrix Factorization: Robust Extraction of Extended Structures
"""
@with_kw struct NMF <: LinearAlgorithm
    ncomps::Union{Int,Nothing} = nothing
end

# manually check for negatives and rescale for full-frame processing
function (alg::NMF)(cube::AbstractArray{T,3}, angles; kwargs...) where T
    target = any(<(0), cube) ? cube .- minimum(cube) : cube 
    S = reconstruct(alg, target, angles)
    residual_cube = target .- S
    return collapse!(residual_cube, angles; kwargs...)
end

# manually check for negatives and rescale for full-frame processing
function (alg::NMF)(cube::AbstractArray{T,3}, angles, cube_ref; kwargs...) where T
    if any(<(0), cube)
        minval = min(minimum(cube), minimum(cube_ref))
        target = cube .- minval
        ref = cube_ref .- minval
    else
        target = cube
        ref = cube_ref
    end
    S = reconstruct(alg, target, angles, ref)
    residual_cube = target .- S
    return collapse!(residual_cube, angles; kwargs...)
end

# TODO this doesn't quite match what sklearn does
function ADI.decompose(alg::NMF, cube, angles, cube_ref; kwargs...)
    if any(<(0), cube) || any(<(0), cube_ref)
        @warn "Negative values encountered in `cube` or `cube_ref`. Make sure to rescale your inputs"
    end
    X = flatten(cube)
    X_ref = flatten(cube_ref)

    k = isnothing(alg.ncomps) ? size(cube, 1) : alg.ncomps
    k > size(cube, 1) && error("ncomps ($k) cannot be greater than the number of frames ($(size(cube, 1)))")
    # fit H using reference
    _, H = nndsvd(X_ref, k)
    W = X * H'
    solve!(CoordinateDescent{eltype(X)}(), X, W, H)
    return H, W
end

function ADI.decompose(alg::NMF, cube, angles; kwargs...)
    # X = clamp.(flatten(cube), 0, Inf) # clamp non-negative values
    if any(<(0), cube)
        @warn "Negative values encountered in `cube`. Make sure to rescale your inputs"
    end
    X = flatten(cube)

    k = isnothing(alg.ncomps) ? size(cube, 1) : alg.ncomps
    k > size(cube, 1) && error("ncomps ($k) cannot be greater than the number of frames ($(size(cube, 1)))")

    W, H = nndsvd(X, k)
    solve!(CoordinateDescent{eltype(X)}(), X, W, H)
    return H, W
end
