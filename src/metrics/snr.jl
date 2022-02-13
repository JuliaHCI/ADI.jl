using Distributions
using HCIToolbox: get_annulus_segments, inside_annulus
using ImageTransformations: center
using Photometry
using Rotations
using StaticArrays

"""
    detectionmap([method=snr], data, fwhm; fill=0)

Parallel implementation of arbitrary detection mapping applied to each pixel in the input image. Any invalid values will be set to `fill`.

The following methods are provided in the [`Metrics`](@ref) module:
* [`snr`](@ref) - signal-to-noise ratio (S/N) using student-t statistics to account for small sample penalty.
* [`significance`](@ref) - Gaussian signifance using student-t statistics to account for small samples penalty.
* [`noise`](@ref) - Standard deviation of apertures in each annulus.

!!! tip
    This code is automatically multi-threaded, so be sure to set `JULIA_NUM_THREADS` before loading your runtime to take advantage of it!
"""
function detectionmap(method, data::AbstractMatrix{T}, fwhm; fill=zero(T)) where T
    out = fill!(similar(data), fill)
    inner_rad = 0.5 * fwhm + 2
    outer_rad = 0.5 * (minimum(size(data)) - fwhm)
    ctr = center(data)

    Threads.@threads for idx in  CartesianIndices(data)
        inside_annulus(inner_rad, outer_rad, ctr, idx) || continue
        val = method(data, idx, fwhm)
        @inbounds out[idx] = ifelse(isfinite(val), val, fill)
    end

    return out
end

detectionmap(data, fwhm) = detectionmap(snr, data, fwhm)

"""
    snr(data, position, fwhm)

Calculate the signal to noise ratio (SNR, S/N) for a test point at `position` using apertures of diameter `fwhm` in a residual frame.

Uses the method of Mawet et al. (2014) which includes penalties for small sample statistics. These are encoded by using a student's t-test for calculating the SNR.

!!! note
    SNR is not equivalent to significance, use [`significance`](@ref) instead
"""
function snr(data::AbstractMatrix, position, fwhm)
    separation = radial_distance(position, center(data))
    separation > fwhm / 2 + 1 || return NaN

    fluxes = get_aperture_fluxes(data, position, separation, fwhm)

    # fluxes of elements which are NOT the test element
    other_elements = @view fluxes[begin + 1:end]
    bkg_μ = mean(other_elements)
    bkg_σ = std(other_elements; corrected=false, mean=bkg_μ)
    sig_factor = sqrt(1 + inv(length(fluxes) - 1))
    return (first(fluxes) - bkg_μ) / (bkg_σ * sig_factor)
end

snr(data::AbstractMatrix, idx::CartesianIndex, fwhm) = snr(data, idx.I, fwhm)

"""
    significance(data, position, fwhm)

Calculates the Gaussian significance from the signal-to-noise ratio (SNR, S/N) for a test point at `position` using apertures of diameter `fwhm` in a residual frame.

The Gaussian signifiance is calculated from converting the SNR confidence limit from a student t distribution to a Gaussian via

``\\text{sig}(\\text{SNR}) = \\Phi^{-1}\\left[\\int_0^\\text{SNR}{t_\\nu(x)dx}\\right]``

where the degrees of freedom ``\\nu`` is given as ``2\\pi r / \\Gamma - 2`` where r is the radial distance of each pixel from the center of the frame.

# See Also
[`snr`](@ref)
"""
function significance(data::AbstractMatrix, position, fwhm)
    separation = radial_distance(position, center(data))
    _snr = snr(data, position, fwhm)
    # put in one line to allow broadcast fusion
    return snr_to_sig(_snr, separation, fwhm)
end
significance(data::AbstractMatrix, idx::CartesianIndex, fwhm) = significance(data, idx.I, fwhm)

function snr_to_sig(snr, separation, fwhm)
    dof = floor(Int, 2 * π * separation / fwhm) - 2
    dof > 0 || return NaN
    @show dof
    return quantile(Normal(), cdf(TDist(dof), float(snr)))
end
function sig_to_snr(sig, separation, fwhm)
    dof = floor(Int, 2 * π * separation / fwhm) - 2
    dof > 0 || return NaN
    return quantile(TDist(dof), cdf(Normal(), float(sig)))
end

"""
    noise(data, position, fwhm)

Calculate the statistical noise for a test point at `position` using apertures of diameter `fwhm` in a residual frame.

Uses the standard deviation of the apertures in the entire annulus. This is distinct from the [`snr`](@ref) noise calculation, which defines a confidence interval using student-t statistics. This means you cannot simply create a noise map and divide it from the signal to create an equivalent S/N map.
"""
function noise(data::AbstractMatrix, position, fwhm)
    separation = radial_distance(position, center(data))
    separation > fwhm / 2 + 1 || return NaN

    fluxes = get_aperture_fluxes(data, position, separation, fwhm)

    return std(fluxes; corrected=false)
end

noise(data::AbstractMatrix, idx::CartesianIndex, fwhm) = noise(data, idx.I, fwhm)

@inline function radial_distance(point, center)
    return sqrt((point[1] - center[1])^2 + (point[2] - center[2])^2)
end

function get_aperture_fluxes(data, position, separation, fwhm)
    r = fwhm / 2

    # number of apertures
    N = floor(Int, 2 * π  * separation / fwhm)
    # rotation per annulus (CW) NB direction doesn't matter?
    dθ = 2 * π / N
    # create affine map rotates points CW by dθ
    transform = recenter(RotMatrix{2}(dθ), center(data))

    # first point
    point = position

    fluxes = similar(data, N)
    @inbounds for idx in eachindex(fluxes)
        # get aperture flux
        ap = CircularAperture(point..., r)
        fluxes[idx] = photometry(ap, data).aperture_sum
        # rotate point
        point = transform(point)
    end

    return fluxes
end
