
using NMF: nndsvd, solve!, CoordinateDescent

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

!!! danger "RDI"
    RDI is not currently supported for NMF due to upstream limitations in NMF.jl

# References
1. [Ren et al. 2018](http://adsabs.harvard.edu/abs/2018ApJ...852..104R) Non-negative Matrix Factorization: Robust Extraction of Extended Structures
"""
struct NMF <: ADIAlgorithm
    ncomps
end
NMF(; ncomps=nothing) = NMF(ncomps)

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

struct NMFDesign{AT,WT} <: LinearDesign
    components::AT
    weights::WT
end

design(des::NMFDesign) = (des.components, des.weights)

function fit(alg::NMF, data::AbstractMatrix{T}) where T
    if any(<(0), data)
        @warn "Negative values encountered in input. Make sure to rescale your data"
    end
    if alg.ncomps === nothing
        k = size(data, 1)
    else
        k = min(alg.ncomps, size(data, 1))
    end

    W, H = nndsvd(data, k)
    solve!(CoordinateDescent{T}(), data, W, H)
    return NMFDesign(H, W)
end

# # TODO this doesn't quite match what sklearn does
# function ADI.decompose(alg::NMF, cube, angles, cube_ref; kwargs...)
#     if any(<(0), cube) || any(<(0), cube_ref)
#         @warn "Negative values encountered in `cube` or `cube_ref`. Make sure to rescale your inputs"
#     end
#     X = flatten(cube)
#     X_ref = flatten(cube_ref)

#     k = isnothing(alg.ncomps) ? size(cube, 1) : alg.ncomps
#     k > size(cube, 1) && error("ncomps ($k) cannot be greater than the number of frames ($(size(cube, 1)))")
#     # fit H using reference
#     _, H = nndsvd(X_ref, k)
#     W = X * H'
#     solve!(CoordinateDescent{eltype(X)}(), X, W, H)
#     return H, W
# end

# function ADI.decompose(alg::NMF, cube, angles; kwargs...)
#     # X = clamp.(flatten(cube), 0, Inf) # clamp non-negative values
#     if any(<(0), cube)
#         @warn "Negative values encountered in `cube`. Make sure to rescale your inputs"
#     end
#     X = flatten(cube)

#     k = isnothing(alg.ncomps) ? size(cube, 1) : alg.ncomps
#     k > size(cube, 1) && error("ncomps ($k) cannot be greater than the number of frames ($(size(cube, 1)))")

#     W, H = nndsvd(X, k)
#     solve!(CoordinateDescent{eltype(X)}(), X, W, H)
#     return H, W
# end
