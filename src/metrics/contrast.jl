using CoordinateTransformations
using Dierckx
using HCIToolbox
using ImageTransformations: center
using LinearAlgebra: dot
using Photometry
using ProgressLogging
using StaticKernels

"""
    contrast_curve(alg, cube, angles, psf;
                   fwhm, sigma=5, nbranch=1, theta=0, inner_rad=1,
                   starphot=Metrics.estimate_starphot(cube, fwhm),
                   fc_rad_sep=3, snr=100, k=2, smooth=true,
                   subsample=true, kwargs...)

Calculate the throughput-calibrated contrast. This first processes the algorithmic [`throughput`](@ref) by injecting instances of `psf` into `cube`. These are processed through `alg` and the ratio of the recovered flux to the injected flux is calculated. These companions are injected in resolution elements across the frame, which can be changed via the various keyword arguments.

The throughput can only be calculated for discrete resolution elements, but we typically want a much smoother curve. To accomplish this, we measure the noise (the standard deviation of all resolution elements in an annulus at a given radius) for every pixel in increasing radii. We then interpolate the throughput to this grid and return the subsampled curves.

# Returned Fields
* `distance` - The radial distance (in pixels) for each measurement
* `contrast` - The Gaussian sensitivity
* `contrast_corr` - The Student-t sensitivity
* `noise` - The noise measured for each distance
* `throughput` - The throughput measured for each distance.

# Keyword Arguments
* `sigma` - The confidence limit in terms of Gaussian standard deviations
* `starphot` - The flux of the star. By default calculates the flux in the central core.
**Injection Options** (See also [`throughput`](@ref))
* `nbranch` - number of azimuthal branches to use
* `theta` - position angle of initial branch
* `inner_rad` - position of innermost planet in FWHM
* `fc_rad_sep` - the separation between planets in FWHM for a single reduction
* `snr` - the target signal to noise ratio of the injected planets
**Subsampling Options** (See also [`Metrics.subsample_contrast`](@ref))
* `subsample` - If true, subsamples the throughput measurements to increase density of curve
* `k` - The order of the BSpline used for subsampling the throughput
* `smooth` - If true, will smooth the subsampled noise measurements with a 2nd order Savitzky-Golay filter

Any additional `kwargs` will be passed to `alg` when it is called.

!!! tip
    If you prefer a tabular format, simply pipe the output of this function into any type supporting the Tables.jl interface, e.g.
    ```
    contrast_curve(alg, cube, angles, psf; fwhm=fwhm) |> DataFrame
    ```
"""
function contrast_curve(alg, cube, angles, psf, args...;
        fwhm, sigma=5, nbranch=1, theta=0, inner_rad=1,
        starphot=Metrics.estimate_starphot(cube, fwhm),
        fc_rad_sep=3, snr=100, k=2, smooth=true,
        subsample=true, kwargs...)
    starphot > 0 || @warn "starphot should be positive but was $starphot"
    # measure the noise and throughput in consecutive resolution elements
    # across azimuthal branches
    @info "Calculating Throughput"
    reduced_empty = alg(cube, angles; kwargs...)

    through, meta = throughput(alg, cube, angles, psf, args...;
                               fwhm=fwhm, inner_rad=inner_rad, fc_rad_sep=fc_rad_sep, theta=theta,
                               nbranch=nbranch, snr=snr, reduced_empty=reduced_empty, kwargs...)

    through_mean = mean(through, dims=2) |> vec

    if subsample
        return subsample_contrast(reduced_empty, meta.distance, through_mean;
                                  fwhm=fwhm, starphot=starphot, sigma=sigma, inner_rad=inner_rad,
                                  theta=theta, smooth=smooth, k=k)
    end

    # calculate common terms once
    unit_contrast = @. meta.noise / (through_mean * starphot)
    # gaussian contrast (invalid values become NaN)
    contrast = calculate_contrast.(sigma, unit_contrast)

    # get correction for small-sample statistics
    sigma_corr = correction_factor.(meta.distance, fwhm, sigma)
    # student-t contrast (invalid values become NaN)
    contrast_corr = calculate_contrast.(sigma_corr, unit_contrast)

    return (distance=meta.distance,
            throughput=through_mean,
            contrast=contrast,
            contrast_corr=contrast_corr,
            noise=meta.noise)
end

function correction_factor(radius, fwhm, sigma)
    n_res_els = 2 * π * radius ÷ fwhm
    ss_corr = sqrt(1 + 1 / (n_res_els - 1))
    return quantile(TDist(n_res_els), cdf(Normal(), sigma)) * ss_corr
end

"""
    Metrics.subsample_contrast(empty_frame, distance, throughput;
                               fwhm, starphot, sigma=5, inner_rad=1,
                               theta=0, smooth=true, k=2)

Helper function to subsample and smooth a contrast curve.

Contrast curves, by definition, are calculated with discrete resolution elements. This can cause contrast curves to have very few points instead of appearing as a continuously measured statistic across the pixels. We alleviate this by sub-sampling the throughput (via BSpline interpolation) across each pixel (instead of each resolution element).

The noise can be found efficiently enough, so rather than interpolate we measure the noise in annuli of width `fwhm` increasing in distance by 1 pixel. We measure this noise in `empty_frame`, which should be a 2D reduced ADI frame.

The noise measurements can be noisy, so a 2nd order Savitzky-Golay filter can be applied via `smooth`. This fits a quadratic polynomial over a window of `fwhm/2` points together to reduce high-frequency jitter.

# Examples

Here is an example which calculates the exact contrast curve in addition to a subsampled curve without re-calculating the throughput.

```julia
cube, angles, psf = # load data

alg = PCA(10)
cc = contrast_curve(alg, cube, angles, psf; fwhm=8.4, subsample=false)
reduced_empty = alg(cube, angles)
sp = Metrics.estimate_starphot(cube, 8.4)
cc_sub = Metrics.subsample_contrast(reduced_empty, cc.distance, cc.throughput;
                                    fwhm=8.4, starphot=sp)
```
"""
function subsample_contrast(empty_frame, distance, throughput;
                            fwhm, starphot, sigma=5, inner_rad=1,
                            theta=0, smooth=true, k=2)
    # measure the noise with high sub-sampling-
    # at every pixel instead of every resolution element
    cx, cy = center(empty_frame)
    radii_subsample = first(distance):1:last(distance) + 1
    sint, cost = sincosd(theta)
    noise_subsample = map(radii_subsample) do r
        x = cx + r * cost
        y = cy + r * sint
        Metrics.noise(empty_frame, (x, y), fwhm)
    end
    through_subsample = Spline1D(distance, throughput, k=k)(radii_subsample)

    if smooth
        window_size = min(length(noise_subsample) - 2, round(Int, 2 * fwhm))
        iseven(window_size) && (window_size += 1)
        width = window_size ÷ 2
        coeffs = savgol_coeffs(window_size, 2) |> Tuple
        smooth_kernel = Kernel{(-width:width,)}(w -> dot(coeffs, Tuple(w)))
        noise_smoothed = map(smooth_kernel, extend(noise_subsample, StaticKernels.ExtensionReplicate()))
    else
        noise_smoothed = noise_subsample
    end

    # calculate common terms once
    unit_contrast = @. noise_smoothed / (through_subsample * starphot)
    # gaussian contrast (invalid values become NaN)
    contrast = calculate_contrast.(sigma, unit_contrast)

    # get correction for small-sample statistics
    sigma_corr = correction_factor.(radii_subsample, fwhm, sigma)
    # student-t contrast (invalid values become NaN)
    contrast_corr = calculate_contrast.(sigma_corr, unit_contrast)

    return (distance=radii_subsample,
            throughput=through_subsample,
            contrast=contrast,
            contrast_corr=contrast_corr,
            noise=noise_smoothed)
end

# simple function to inline the contrast calculation _and_ clipping
@inline function calculate_contrast(k, unit_contrast)
    contrast = k * unit_contrast
    return ifelse(0 ≤ contrast ≤ 1, contrast, oftype(contrast, NaN))
end

function savgol_coeffs(window_size, order=2)
    pos = window_size ÷ 2
    # form vandermonde matrix
    x = -pos:window_size - pos - 1
    powers = 0:order
    A = x' .^ powers
    # solve least squares equation Ac = y
    y = zeros(order + 1)
    y[1] = 1
    return A \ y
end

"""
    Metrics.estimate_starphot(cube, fwhm)
    Metrics.estimate_starphot(frame, fwhm)

Simple utility to estimate the stellar photometry by placing a circular aperture with `fwhm` diameter in the center of the `frame`. If a cube is provided, first the median frame will be found.
"""
function estimate_starphot(frame::AbstractMatrix, fwhm)
    ap = CircularAperture(center(frame), fwhm / 2)
    return photometry(ap, frame).aperture_sum
end

estimate_starphot(cube::AbstractArray{T,3}, fwhm) where {T} = estimate_starphot(collapse(cube, method=median), fwhm)

"""
    throughput(alg, cube, angles, psf;
               fwhm, nbranch=1, theta=0, inner_rad=1,
               fc_rad_sep=3, snr=100, kwargs...)

Calculate the throughput of `alg` by injecting fake companions into `cube` and measuring the relative photometry of each companion in the reduced frame. The photometry is measured using a circular aperture with a diameter matching the `fwhm`. Any additional `kwargs` will be passed to `alg` when it is called.

# Keyword Arguments
* `nbranch` - number of azimuthal branches to use
* `theta` - position angle of initial branch
* `inner_rad` - position of innermost planet in FWHM
* `fc_rad_sep` - the separation between planets in FWHM for a single reduction
* `snr` - the target signal to noise ratio of the injected planets
* `reduced_empty` - the collapsed residual frame for estimating the noise. Will process using `alg` if not provided.
"""
function throughput(alg, cube::AbstractArray{T,3}, angles, psf_model;
                    fwhm, nbranch=1, theta=0, inner_rad=1, fc_rad_sep=3,
                    snr=100, reduced_empty = nothing, kwargs...) where T
    maxfcsep = minimum(Base.front(size(cube))) ÷ (2 * fwhm) - 1
    # too large separation between companions in the radial patterns
    3 ≤ fc_rad_sep ≤ maxfcsep || error("`fc_rad_sep` should lie ∈[3, $(maxfcsep)], got $fc_rad_sep")

    # compute noise in concentric annuli on the empty frame
    if isnothing(reduced_empty)
        reduced_empty = alg(cube, angles; kwargs...)
    end

    cx, cy = center(reduced_empty)

    n_annuli = minimum(floor.(Int, ((cx, cy) .- fwhm) ./ fwhm) .- 1)
    radii = _fix_range(fwhm .* (inner_rad:n_annuli), cube)

    sint, cost = sincosd(theta)
    noise = map(radii) do r
        x = cx + r * cost
        y = cy + r * sint
        Metrics.noise(reduced_empty, (x, y), fwhm)
    end

    angle_per_branch = 360 / nbranch
    output = similar(cube, length(radii), nbranch)

    tmp_cube = similar(cube)
    fake_comps_full = zero(reduced_empty)
    tmp_frame = similar(reduced_empty)
    @progress name="branch" for branch in 1:nbranch
        θ = theta + angle_per_branch * (branch - 1)
        δy, δx = sincosd(θ)
        @progress name="pattern" for init_rad in 1:fc_rad_sep
            slice = init_rad:fc_rad_sep:lastindex(radii)
            
            # reset work arrays
            fill!(tmp_frame, zero(T))
            copyto!(tmp_cube, cube)

            apertures = map(slice) do ann
                r = radii[ann]
                x = r * δx + cx
                y = r * δy + cy

                A = snr * noise[ann]

                inject!(tmp_frame, psf_model; x, y, amp=A, fwhm)
                inject!(tmp_cube, psf_model, angles; x, y, amp=A, fwhm)

                return CircularAperture(x, y, fwhm / 2)
            end
            # add fake comps to output frame
            fake_comps_full .+= tmp_frame

            # get reduced output after subtracting and collapsing
            reduced = alg(tmp_cube, angles; kwargs...)

            injected_flux = photometry(apertures, tmp_frame).aperture_sum
            recovered_flux = photometry(apertures, reduced .- reduced_empty).aperture_sum
            @. output[slice, branch] = max(zero(T), recovered_flux / injected_flux)
        end
    end

    return output, (distance=radii, fake_comps=fake_comps_full, noise=noise)
end

_fix_range(radii, cube::AbstractArray{T,3}) where {T} = radii
_fix_range(radii, cube::AnnulusView) = filter(r -> cube.rmin ≤ r ≤ cube.rmax, radii)
function _fix_range(radii, cube::MultiAnnulusView)
    rmin, rmax = extrema(cube.radii) .- cube.width / 2
    filter(r -> rmin ≤ r ≤ rmax, radii)
end

"""
    throughput(alg, cube, angles, psf, position;
               fwhm, snr=100,
               reduced_empty=nothing,
               verbose=true, kwargs...)

Calculate the throughput of `alg` by injecting `psf` into `cube` at the given `position` and measuring the relative photometry of the companion in the reduced frame. The photometry is measured using a circular aperture with a diameter matching the `fwhm`. Any additional `kwargs` will be passed to `alg` when it is called.

If `position` is a tuple or a vector, it will be parsed as the cartesian coordinate `(x, y)`. If `position` is a `CoordinateTransformations.Polar` it will be parsed as polar coordinates from the center of the cube. Note the `Polar` type expects angles in radians.

# Keyword Arguments
* `fwhm` - the full width at half-maximum
* `snr` - the target signal to noise ratio of the injected planets
* `reduced_empty` - the collapsed residual frame for estimating the noise. Will process using `alg` if not provided.
* `verbose` - show informative messages during process
"""
function throughput(alg, cube::AbstractArray{T,3}, angles, psf_model, (x, y); fwhm,
    snr=100, reduced_empty=nothing, verbose=true, kwargs...) where T

    if reduced_empty === nothing
        verbose && @info "Calculating empty reduced frame"
        reduced_empty = alg(cube, angles; kwargs...)
    end

    # find amplitude of injected PSF
    noise = Metrics.noise(reduced_empty, (x, y), fwhm)
    A = snr * noise

    verbose && @info "Injecting companion at x=$x y=$y with amp=$A"
    fake_comp = inject!(zero(reduced_empty), psf_model; x, y, amp=A, fwhm)
    fake_comp_cube = inject(cube, psf_model, angles; x, y, amp=A, fwhm)

    verbose && @info "Calculating reduced frame with fake companion injected"
    reduced = alg(fake_comp_cube, angles; kwargs...)

    ap = CircularAperture(x, y, fwhm/2)

    verbose && @info "Calculating throughput"
    injected_flux = photometry(ap, fake_comp).aperture_sum
    recovered_flux = photometry(ap, reduced .- reduced_empty).aperture_sum
    throughput = max(zero(T), recovered_flux / injected_flux)

    return throughput, (;fake_comp, fake_comp_cube, reduced, injected_flux, recovered_flux)
end
