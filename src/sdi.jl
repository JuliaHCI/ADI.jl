using FillArrays

"""
    ADI.SDIAlgorithm <: ADI.ADIAlgorithm

Spectral differential imaging (SDI) algorithms. These work on 4-D SDI tensors. To use these algorithms, simply treat them like functions

```julia
(::SDIAlgorithm)(data::AbstractArray{T,4}, angles, scales; [ref] kwargs...)
```

The data is expected to be laid out in `(nλ, nf, ny, nx)` format, so you may need to `permutedims` before processing the data. The `scales` correspond to the relative wavelength scales for each spectrum, and can be retrieved with `HCIToolbox.scale_list`.

# Algorithms

The current SDI implementations are
* [`SingleSDI`](@ref)
* [`DoubleSDI`](@ref)
* [`SliceSDI`](@ref)
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

function (sdi::SingleSDI)(spcube::AbstractArray{T,4}, angles, scales; method=:deweight, kwargs...) where T
    nλ, n, ny, nx = size(spcube)
    frame_size = (ny, nx)
    big_cube = scale_and_stack(spcube, scales)
    angs = repeat(angles, inner=nλ)
    # do single-pass reconstruction
    if :ref in keys(kwargs)
        big_cube_ref = scale_and_stack(kwargs[:ref], scales)
        big_resid_cube = subtract(sdi.alg, big_cube; angles=angs, kwargs..., ref=big_cube_ref)
    else
        big_resid_cube = subtract(sdi.alg, big_cube; angles=angs, kwargs...)
    end

    # bin across spectral dim
    resid_cube = invscale_and_collapse(big_resid_cube, scales, frame_size)
    # derotate and combine
    return collapse!(resid_cube, angles; method=method, kwargs...)
end

"""
    DoubleSDI(alg)
    DoubleSDI(alg_spec, alg_temp)

A wrapper algorithm for spectral differential imaging (SDI) data reduced in two passes. The first pass uses `alg_spec` to reduce each spectral cube slice in the SDI tensor. Then, the spectral residual frames will be reduced using `alg_temp`, which will include the derotation and final combination.

The difference between [`DoubleSDI`](@ref) and [`SliceSDI`](@ref) is that `DoubleSDI` does its first pass in the spectral slice, effectively collapsing the slice before performing ADI on the residual cube. `SliceSDI` does its first pass in the temporal slice, collapsing it first before performing ADI on the residual cube.

# Examples

```julia
julia> data, angles, scales = # load data...

# Median subtraction for each spectral slice,
# GreeDS{PCA} subtraction on spectral residual cube
julia> res = DoubleSDI(Classic(), GreeDS(15))(data, angles, scales)
```
"""
struct DoubleSDI{ALG<:ADIAlgorithm, ALG2<:ADIAlgorithm} <: SDIAlgorithm
    alg_spec::ALG
    alg_temp::ALG2
end

DoubleSDI(alg) = DoubleSDI(alg, alg)

function (sdi::DoubleSDI)(spcube::AbstractArray{T,4}, angles, scales; method=:deweight, kwargs...) where T
    nλ, n, ny, nx = size(spcube)
    frame_size = (ny, nx)
    spec_resids = similar(spcube, n, ny, nx)
    # do first pass in spectral domain
    Threads.@threads for n in axes(spcube, 2)
        cube = @view spcube[:, n, :, :]
        scaled_cube = scale(cube, scales)
        angs = Zeros(nλ)
        if :ref in keys(kwargs)
            cube_ref = @view kwargs[:ref][:, n, :, :]
            scaled_cube_ref = scale(cube_ref, scales)
            R = subtract(sdi.alg_spec, scaled_cube; angles=angs, kwargs..., ref=scaled_cube_ref)
        else
            R = subtract(sdi.alg_spec, scaled_cube; angles=angs, kwargs...)
        end
        spec_resid = invscale(R, scales, frame_size)
        spec_resids[n, :, :] .= collapse(spec_resid)
    end
    # do second pass in temporal domain
    return sdi.alg_temp(spec_resids, angles; method=method, angles=angles, kwargs...)
end

"""
    SliceSDI(alg)
    SliceSDI(alg_spec, alg_temp)

A wrapper algorithm for spectral differential imaging (SDI) data reduced in two passes. The first pass uses `alg_temp` to reduce each temporal cube slice in the SDI tensor. These residuals will be rescaled and stacked into a new cube. Then, the temporal residual frames will be reduced using `alg_spec`, which will include the derotation and final combination.

The difference between [`SliceSDI`](@ref) and [`DoubleSDI`](@ref) is that `DoubleSDI` does its first pass in the spectral slice, effectively collapsing the slice before performing ADI on the residual cube. `SliceSDI` does its first pass in the temporal slice, collapsing it first before performing ADI on the residual cube.

# Examples

```julia
julia> data, angles, scales = # load data...

# Median subtraction for each spectral slice,
# GreeDS{PCA} subtraction on spectral residual cube
julia> res = SliceSDI(Classic(), GreeDS(15))(data, angles, scales)
```
"""
struct SliceSDI{ALG<:ADIAlgorithm, ALG2<:ADIAlgorithm} <: SDIAlgorithm
    alg_spec::ALG
    alg_temp::ALG2
end

SliceSDI(alg) = SliceSDI(alg, alg)

function (sdi::SliceSDI)(spcube::AbstractArray{T,4}, angles, scales; kwargs...) where T
    nλ, n, ny, nx = size(spcube)
    frame_size = (ny, nx)
    temp_resids = similar(spcube, nλ, ny, nx)
    # do first pass in temporal domain
    Threads.@threads for n in axes(spcube, 1)
        cube = spcube[n, :, :, :]
        if :ref in keys(kwargs)
            cube_ref = kwargs[:ref][n, :, :, :]
            temp_resids[n, :, :] .= sdi.alg_temp(cube, angles; kwargs..., ref=cube_ref)
        else
            temp_resids[n, :, :] .= sdi.alg_temp(cube, angles; kwargs...)
        end
    end
    # do second pass in temporal domain
    scaled_resid_cube = scale(temp_resids, scales)
    angs = Zeros(nλ)
    resid = sdi.alg_spec(scaled_resid_cube, angs; angles=angs, kwargs...)
    return invscale(resid, maximum(scales))
end
