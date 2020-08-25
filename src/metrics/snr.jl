using ImageTransformations: center
using Photometry
using Distributions
using Statistics
using ImageFiltering
using StatsBase: mad
using HCIToolbox: get_annulus_segments

"""
    snrmap(data, fwhm; method=:exact)

Parallel implementation of signal-to-noise ratio (SNR, S/N) applied to each pixel in the input image.



### SNR Methods
* `:exact` - Uses [`snr`](@ref) to get the exact SNR using student t statistics (small samples penalty)
* `:approx` - Convolves the input by a tophat kernel before using [`snr_approx!`](@ref) to approximate the student t statistics (small samples penalty).

!!! tip
    This code is automatically multi-threaded, so be sure to set `JULIA_NUM_THREADS` before loading your runtime to take advantage of it!
"""
function snrmap(data::AbstractMatrix{T}, fwhm; method=:exact) where T
    out = fill!(similar(data), zero(T))
    width = minimum(size(data)) / 2 - 1.5 * fwhm

    data = _prepmatrix(Val(method), data, fwhm)

    masked = get_annulus_segments(data, fwhm/2 + 2, width, mode=:apply)
    coords = findall(!iszero, masked)

    Threads.@threads for coord in coords
        @inbounds out[coord] = _snrfunc(Val(method))(data, coord, fwhm)
    end
    
    return out
end

"""
    snr(data, position, fwhm)

Calculate the signal to noise ratio (SNR, S/N) for a test point at `position` using apertures of diameter `fwhm` in a residual frame.

Uses the method of Mawet et al. 2014 which includes penalties for small sample statistics. These are encoded by using a student's t-test for calculating the SNR.

!!! note
    SNR is not equivalent to significance, use [`significance`](@ref) instead
"""
function snr(data::AbstractMatrix, position, fwhm)
    x, y = position
    cy, cx = center(data)
    separation = sqrt((x - cx)^2 + (y - cy)^2)
    @assert separation > fwhm / 2 + 1 "`position` is too close to the frame center"

    θ = 2asin(fwhm / 2 / separation)
    N = floor(Int, 2π / θ)

    sint, cost = sincos(θ)
    xs = similar(data, N)
    ys = similar(data, N)

    # initial points
    rx = x - cx
    ry = y - cy

    @inbounds for idx in eachindex(xs)
        xs[idx] = rx + cx
        ys[idx] = ry + cy
        rx, ry = cost * rx + sint * ry, cost * ry - sint * rx
    end

    r = fwhm / 2

    apertures = CircularAperture.(xs, ys, r)
    fluxes = aperture_photometry(apertures, data, method=:exact).aperture_sum
    other_elements = @view fluxes[2:end]
    bkg_σ = std(other_elements) # ddof = 1 by default
    return (fluxes[1] - mean(other_elements)) / (bkg_σ * sqrt(1 + 1/(N - 1)))
end

snr(data::AbstractMatrix, idx::CartesianIndex, fwhm) = snr(data, (idx.I[2], idx.I[1]), fwhm)


"""
    snr_approx!(data, position, fwhm)

In-place version of [`snr_approx`](@ref) which modifies `data`.
"""
function snr_approx!(data::AbstractMatrix, position, fwhm)
    x, y = position
    cy, cx = center(data)
    separation = sqrt((x - cx)^2 + (y - cy)^2)
    
    aper_ind = circle_index((y, x), fwhm/2)
    ann_ind = annulus_index(floor.(Int, (cy, cx)), floor(Int, separation))

    peak = data[y, x]

    data[aper_ind] .= mad(data[ann_ind], normalize=false)
    N = 2π * separation / fwhm - 1
    noise = std(data[ann_ind], corrected=false) * sqrt(1 + 1/N)
    signal = peak - mean(data[ann_ind])
    return signal / noise
end

snr_approx!(data::AbstractMatrix, idx::CartesianIndex, fwhm) = snr_approx!(data, (idx.I[2], idx.I[1]), fwhm)

"""
    snr_approx(data, position, fwhm)

Calculates an approximate signal to noise ratio (SNR, S/N) for a test point at `position` using apertures of diameter `fwhm` in a residual frame.

Applies the small sample statistics penalty the same as [`snr`](@ref). Data is assumed to have been filtered using a 2D top-hat kernel already (automatically done if called via [`snrmap`](@ref))

!!! note
    SNR is not equivalent to significance, use [`significance`](@ref) instead.
"""
snr_approx(data, position, fwhm) = snr_approx!(deepcopy(data), position, fwhm)

## These methods need to be below function definitions
# no prep needed for exact method
_prepmatrix(::Val{:exact}, data, fwhm) = data
# applies a 2D top-hat filter with radius fwhm/2
# TODO watch for imagefiltering merge and switch when that happens
function _prepmatrix(::Val{:approx}, data, fwhm)
    sz = _round_up_to_odd_integer(fwhm)
    kern = similar(data, sz, sz)
    ctr = center(kern)
    for idx in CartesianIndices(kern)
        d = sqrt((idx[1] - ctr[1])^2 + (idx[2] - ctr[2])^2)
        kern[idx] = d < fwhm/2 ? 1 : 0
    end
    kern ./= sum(kern)
    return imfilter(data, centered(kern), Fill(0))
end

function _round_up_to_odd_integer(value)
    i = ceil(Int, value)
    return iseven(i) ? i + 1 : i
end


_snrfunc(::Val{:exact}) = snr
_snrfunc(::Val{:approx}) = snr_approx!

## 


"""
    significance(snr::AbstractMatrix, fwhm)

Calculates the Gaussian significance from the signal-to-noise ratio (SNR, S/N).

The Gaussian signifiance is calculated from converting the SNR confidence limit from a student t distribution to a Gaussian via

``\\text{sig}(\\text{SNR}) = \\Phi^{-1}\\left[\\int_0^\\text{SNR}{t_\\nu(x)dx}\\right]``

where the degrees of freedom ``\\nu`` is given as ``2\\pi r / \\Gamma - 2`` where r is the radial distance of each pixel from the center of the frame.

# See Also
[`snrmap`](@ref)
"""
function significance(snr::AbstractMatrix, fwhm)
    ys, xs = axes(snr)
    cy, cx = center(snr)
    # put in one line to allow broadcast fusion
    return @. snr_to_sig(snr, sqrt((xs' - cx)^2 + (ys - cy)^2), fwhm)
end

snr_to_sig(snr, separation, fwhm) = @. quantile(Normal(), cdf(TDist(2π * separation / fwhm .- 2), Float64(snr)))
sig_to_snr(sig, separation, fwhm) = @. quantile(TDist(2π * separation / fwhm.- 2), cdf(Normal(), Float64(sig)))

# equivalent to skimage.draw.circle
function circle_index(center, r)
    upper_left = @. ceil(Int, center - r)
    lower_right = @. floor(Int, center + r)
    shape = @. lower_right - upper_left
    shift_center = @. center - upper_left
    rows = Base.OneTo(shape[1]) .- shift_center[1]
    cols = Base.OneTo(shape[2]) .- shift_center[2]
    dist = @. (rows / r)^2 + (cols' / r)^2
    idxs = findall(v -> v < 1, dist)
    return [idx + CartesianIndex(upper_left...) for idx in idxs]
end

# equivalent to skimage.draw.circle_perimiter using the Bresenham method
function annulus_index(center, r)
    d = 3 - 2r
    
    row, col = r, 0
    rows = Int[]
    cols = Int[]
    while row ≥ col
        push!(rows, row, -row, row, -row, col, -col, col, -col)
        push!(cols, col, col, -col, -col, row, row, -row, -row)

        if d < 0
            d += 4col + 6
        else
            d += 4 * (col - row) + 10
            row -= 1
        end
        col += 1
    end
    return @. CartesianIndex(rows + center[1], cols + center[2])
end
