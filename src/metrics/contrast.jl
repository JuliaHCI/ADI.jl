using Statistics
using ImageTransformations: center
using Photometry
using HCIToolbox
using Dierckx
using ProgressLogging

"""
    contrast_curve(alg,
                   cube,
                   angles,
                   psf,
                   args...;
                   fwhm,
                   sigma=5,
                   starphot=Metrics.estimate_starphot(cube, fwhm),
                   nbranch=1,
                   theta=0,
                   inner_rad=1,
                   fc_rad_sep=3,
                   snr=100,
                   k=2,
                   kwargs...)

Calculate the throughput-calibrated contrast. This first processes the algorithmic [`throughput`](@ref) by injecting instances of `psf` into `cube`. These are processed through `alg` and the ratio of the recovered flux to the injected flux is calculated. These companions are injected in resolution elements across the frame, which can be changed via the various keyword arguments.

The throughput can only be calculated for discrete resolution elements, but we typically want a much smoother curve. To accomplish this, we measure the noise (the standard deviation of all resolution elements in an annulus at a given radius) for every pixel in increasing radii. We then interpolate the throughput to this grid and return the subsampled curves.

# Fields
* `distance` - The radial distance (in pixels) for each measurement
* `contrast` - The Gaussian sensitivity
* `contrast_corr` - The Student-t sensitivity
* `noise` - The noise measured for each distance
* `throughput` - The throughput measured for each distance.

# Keyword Arguments
* `sigma` - The confidence limit in terms of Gaussian standard deviations
* `starphot` - The flux of the star. By default calculates the flux in the central core.
* `nbranch` - number of azimuthal branches to use
* `theta` - position angle of initial branch
* `inner_rad` - position of innermost planet in FWHM
* `fc_rad_sep` - the separation between planets in FWHM for a single reduction
* `snr` - the target signal to noise ratio of the injected planets
* `k` - The order of the BSpline used for subsampling the throughput

!!! tip
    If you prefer a tabular format, simply pipe the output of this function into any type supporting the Tables.jl interface, e.g.
    ```
    contrast_curve(alg, cube, angles, psf; fwhm=fwhm) |> DataFrame
    ```
"""
function contrast_curve(alg,
                        cube,
                        angles,
                        psf,
                        args...;
                        fwhm,
                        starphot=estimate_starphot(cube, fwhm),
                        sigma=5,
                        inner_rad=1,
                        k=2,
                        kwargs...)

    # measure the noise and throughput in consecutive resolution elements
    # across azimuthal branches
    @info "Calculating Throughput"
    through, meta = throughput(alg, cube, angles, psf, args...; fwhm=fwhm, inner_rad=inner_rad, kwargs...)

    through_mean = mean(through, dims=2) |> vec

    @info "Calculating pixel-wise noise"
    # measure the noise with high sub-sampling, at every pixel instead of every resolution element
    reduced_empty = alg(cube, angles, args...)
    _, cy, cx = center(cube)
    radii_subsample = (fwhm * inner_rad):(cy - fwhm / 2)
    noise_subsample = map(radii_subsample) do r
        x = r + cx
        y = cy
        Metrics.noise(reduced_empty, (x, y), fwhm)
    end

    @info "Calculating contrast"

    through_subsample = Spline1D(meta.distance, through_mean, k=k)(radii_subsample)

    unit_contrast = @. noise_subsample / (through_subsample * starphot)
    contrast = @. sigma * unit_contrast
    @. contrast[!(0 ≤ contrast ≤ 1)] = NaN

    # get correction for small-sample statistics
    sigma_corr = @. correction_factor.(radii_subsample, fwhm, sigma)

    contrast_corr = @. sigma_corr * unit_contrast
    @. contrast_corr[!(0 ≤ contrast_corr ≤ 1)] = NaN

    return (distance=radii_subsample,
            throughput=through_subsample,
            contrast=contrast,
            contrast_corr=contrast_corr,
            noise=noise_subsample)
end

function correction_factor(radius, fwhm, sigma)
    n_res_els = 2 * π * radius ÷ fwhm
    ss_corr = sqrt(1 + 1 / (n_res_els - 1))
    return quantile(TDist(n_res_els), cdf(Normal(), sigma)) * ss_corr
end

"""
    Metrics.estimate_starphot(cube, fwhm)
    Metrics.estimate_starphot(frame, fwhm)

Simple utility to estimate the stellar photometry by placing a circular aperture with `fwhm` diameter in the center of the `frame`. If a cube is provided, first the median frame will be found.
"""
function estimate_starphot(frame::AbstractMatrix, fwhm)
    ap = CircularAperture(reverse(center(frame)), fwhm/2)
    return photometry(ap, frame).aperture_sum
end

estimate_starphot(cube::AbstractArray{T, 3}, fwhm) where {T} = estimate_starphot(collapse(cube, method=median), fwhm)


"""
    throughput(alg,
               cube,
               angles,
               psf,
               args...;
               fwhm,
               nbranch=1,
               theta=0,
               inner_rad=1,
               fc_rad_sep=3,
               snr=100,
               kwargs...)

Calculate the throughput of `alg` by injecting fake companions into `cube` and measuring the relative photometry of each companion in the reduced frame. Any additional `args` or `kwargs` will be passed to `alg` when it is called.

# Keyword Arguments
* `nbranch` - number of azimuthal branches to use
* `theta` - position angle of initial branch
* `inner_rad` - position of innermost planet in FWHM
* `fc_rad_sep` - the separation between planets in FWHM for a single reduction
* `snr` - the target signal to noise ratio of the injected planets
"""
function throughput(alg,
                    cube::AbstractArray{T,3},
                    angles,
                    psf_model,
                    args...;
                    fwhm,
                    nbranch=1,
                    theta=0,
                    inner_rad=1,
                    fc_rad_sep=3,
                    snr=100,
                    kwargs...) where T
    maxfcsep = size(cube, 2) ÷ (2 * fwhm) - 1
    # too large separation between companions in the radial patterns
    3 ≤ fc_rad_sep ≤ maxfcsep || error("`fc_rad_sep` should lie ∈[3, $(maxfcsep)], got $fc_rad_sep")

    # compute noise in concentric annuli on the empty frame
    reduced_empty = alg(cube, angles, args...; kwargs...)

    cy, cx = center(reduced_empty)

    n_annuli = floor(Int, (cy - fwhm) / fwhm) - 1
    radii = fwhm .* (inner_rad:n_annuli)
    δy, δx = sincosd(theta)
    noise = map(radii) do r
        x = r * δx + cx
        y = r * δy + cy
        Metrics.noise(reduced_empty, (x, y), fwhm)
    end

    angle_per_branch = 360 / nbranch
    output = similar(cube, length(radii), nbranch)

    fake_comps_full = zero(reduced_empty)
    @progress "branch" for branch in 1:nbranch
        θ = theta + angle_per_branch * (branch - 1)
        @progress "radial slice" for init_rad in 1:fc_rad_sep
            slice = init_rad:fc_rad_sep:lastindex(radii)
            fake_comps = zero(reduced_empty)

            cube_fake_comps = copy(cube)

            apertures = map(slice) do ann
                r = radii[ann]
                δy, δx = sincosd(θ)
                x = r * δx + cx
                y = r * δy + cy

                A = snr * noise[ann]

                inject!(fake_comps, psf_model; A=A, r=r, θ=θ)
                fake_comps_full .+= fake_comps
                inject!(cube_fake_comps, psf_model, angles; A=A, r=r, θ=θ)

                return CircularAperture(x, y, fwhm / 2)
            end
            reduced = alg(cube_fake_comps, angles, args...; kwargs...)

            injected_flux = photometry(apertures, fake_comps).aperture_sum
            recovered_flux = photometry(apertures, reduced .- reduced_empty).aperture_sum
            @. output[slice, branch] = max(zero(T), recovered_flux / injected_flux)
        end
    end

    return output, (distance=radii, fake_comps=fake_comps_full, noise=noise)
end
