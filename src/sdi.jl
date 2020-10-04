using FillArrays

"""
    ADI.SDIAlgorithm <: ADI.ADIAlgorithm

Spectral differential imaging (SDI) algorithms. These work on 4-D SDI tensors. To use these algorithms, simply treat them like functions

```julia
(::SDIAlgorithm)(data::AbstractArray{T,4}, angles, scales; kwargs...)
(::SDIAlgorithm)(data::AbstractArray{T,4}, angles, data_ref, scales; kwargs...)
```

# Algorithms

The current SDI implementations are
* [`SingleSDI`](@ref)
* [`DoubleSDI`](@ref)
"""
abstract type SDIAlgorithm <: ADIAlgorithm end

"""
    SingleSDI(alg)

A wrapper algorithm for spectral differential imaging (SDI) data reduced in a single pass. This means that each channel will be scaled and then concatenated together, so an SDI tensor `(nλ, nf, y, x)` becomes a stack `(nλ * nf, y, x)` which is processed by the underlying `alg` like ADI data.

!!! tip
    `SingleSDI` is the default SDI mode. This means instead of writing
    ```julia
    SingleSDI(PCA(15))(data, angles, scales)
    ```
    you can just write
    ```julia
    PCA(15)(data, angles, scales)
    ```
"""
struct SingleSDI{ALG<:ADIAlgorithm} <: SDIAlgorithm
    alg::ALG
end

(alg::ADIAlgorithm)(spcube::AbstractArray{T,4}, angles, scales; kwargs...) where {T} =
    SingleSDI(alg)(spcube, angles, scales; kwargs...)

(alg::ADIAlgorithm)(spcube::AbstractArray{T,4}, angles, scales, spcube_ref; kwargs...) where {T} =
    SingleSDI(alg)(spcube, angles, scales, spcube_ref; kwargs...)

function (sdi::SingleSDI)(spcube::AbstractArray{T,4}, angles, scales; kwargs...) where T
    nλ, n, ny, nx = size(spcube)
    frame_size = (ny, nx)
    big_cube = scale_and_stack(spcube, scales)

    # do single-pass reconstruction
    S = reconstruct(sdi.alg, big_cube, repeat(angles, inner=nλ); kwargs...)
    big_resid_cube = big_cube .- S
    # bin across spectral dim
    resid_cube = invscale_and_collapse(big_resid_cube, scales, frame_size)
    # derotate and combine
    return collapse!(resid_cube, angles; kwargs...)
end

function (sdi::SingleSDI)(spcube::AbstractArray{T,4}, angles, scales, spcube_ref; kwargs...) where T
    nλ, n, ny, nx = size(spcube)
    frame_size = (ny, nx)
    big_cube = scale_and_stack(spcube, scales)
    big_cube_ref = scale_and_stack(spcube_ref, scales)

    # do single-pass reconstruction
    S = reconstruct(sdi.alg, big_cube, repeat(angles, inner=nλ), big_cube_ref; kwargs...)
    big_resid_cube = big_cube .- S
    # bin across spectral dim
    resid_cube = invscale_and_collapse(big_resid_cube, scales, frame_size)
    # derotate and combine
    return collapse!(resid_cube, angles; kwargs...)
end

"""
    DoubleSDI(alg)
    DoubleSDI(alg_spec, alg_temp)

A wrapper algorithm for spectral differential imaging (SDI) data reduced in two passes. The first pass uses `alg_spec` to reduce each spectral cube slice in the SDI tensor. Then, the spectral residual frames will be reduced using `alg_temp`, which will include the derotation and final combination.

# Examples

```julia
julia> data, angles, scales = # load data...

# Median subtraction for each spectral slice,
# GreeDS{PCA} subtraction on spectral residual cube
julia> res = DoubleSDI(Median(), GreeDS(15))(data, angles, scales)
```
"""
struct DoubleSDI{ALG<:ADIAlgorithm, ALG2<:ADIAlgorithm} <: SDIAlgorithm
    alg_spec::ALG
    alg_temp::ALG2
end

DoubleSDI(alg) = DoubleSDI(alg, alg)

function (sdi::DoubleSDI)(spcube::AbstractArray{T,4}, angles, scales; kwargs...) where T
    nλ, n, ny, nx = size(spcube)
    frame_size = (ny, nx)
    spec_resids = similar(spcube, n, ny, nx)
    # do first pass in spectral domain
    Threads.@threads for n in axes(spcube, 2)
        cube = @view spcube[:, n, :, :]
        scaled_cube = scale(cube, scales)
        S = reconstruct(sdi.alg_spec, scaled_cube, Fill(angles[n], nλ); kwargs...)
        spec_resid = invscale(scaled_cube .- S, scales, frame_size)
        spec_resids[n, :, :] .= collapse(spec_resid)
    end
    # do second pass in temporal domain
    return sdi.alg_temp(spec_resids, angles; kwargs...)
end

function (sdi::DoubleSDI)(spcube::AbstractArray{T,4}, angles, scales, spcube_ref; kwargs...) where T
    nλ, n, ny, nx = size(spcube)
    frame_size = (ny, nx)
    spec_resids = similar(spcube, n, ny, nx)
    # do first pass in spectral domain
    Threads.@threads for n in axes(spcube, 2)
        cube = @view spcube[:, n, :, :]
        cube_ref = @view spcube_ref[:, n, :, :]
        scaled_cube = scale(cube, scales)
        scaled_cube_ref = scale(cube_ref, scales)
        S = reconstruct(sdi.alg_spec, scaled_cube, Fill(angles[n], nλ), scaled_cube_ref; kwargs...)
        spec_resid = invscale(scaled_cube .- S, scales, frame_size)
        spec_resids[n, :, :] .= collapse(spec_resid)
    end
    # do second pass in temporal domain
    return sdi.alg_temp(spec_resids, angles; kwargs...)
end