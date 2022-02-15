
using NMF: nndsvd, solve!, CoordinateDescent

"""
    NMF(;ncomps=nothing)
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

# References
1. [Ren et al. 2018](http://adsabs.harvard.edu/abs/2018ApJ...852..104R) Non-negative Matrix Factorization: Robust Extraction of Extended Structures
"""
@concrete struct NMF <: ADIAlgorithm
    ncomps
end
NMF(; ncomps=nothing) = NMF(ncomps)

function subtract(alg::NMF, cube; kwargs...)
    target = normalize_nmf_input(cube)
    if :ref in keys(kwargs)
        ref_ = normalize_nmf_input(kwargs[:ref])
        S = reconstruct(alg, target; kwargs..., ref=ref_)
    else
        S = reconstruct(alg, target; kwargs...)
    end
    return target .- S
end

function normalize_nmf_input(data)
    if any(<(0), data)
        return data .- minimum(data)
    end
    return data
end

function fit(alg::NMF, data::AbstractMatrix{T}; ref=data, kwargs...) where T
    if any(<(0), data) || any(<(0), ref)
        @warn "Negative values encountered in input. Make sure to rescale your data"
    end
    if alg.ncomps === nothing
        k = size(data, 2)
    else
        k = min(alg.ncomps, size(data, 2))
    end
    # transpose data to accommadate NMF.jl interface
    X = collect(transpose(ref))
    W, H = nndsvd(X, k)
    solve!(CoordinateDescent{T}(), X, W, H)
    if ref !== data
        Y = collect(transpose(data))
        W = Y * transpose(H)
        # have to transpose to use `update_H=false`
        solve!(CoordinateDescent{T}(update_H=false), Y, W, H)
    end
    return LinearDesign(transpose(H), transpose(W))
end
