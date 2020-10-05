using Statistics
using ImageTransformations: center
using Photometry
using HCIToolbox

"""
    contrast_curve
"""
function contrast_curve(alg,
                        cube,
                        angles,
                        psf,
                        args...;
                        fwhm,
                        starphot=1,
                        sigma=5,
                        nbranch=1,
                        theta=0,
                        inner_rad=1,
                        fc_snr=100,
                        transmission=nothing)

    nothing
end

"""
    throughput
"""
function throughput(alg,
                    cube,
                    angles,
                    psf_model;
                    fwhm,
                    sigma=5,
                    nbranch=1,
                    theta=0,
                    inner_rad=1,
                    fc_rad_sep=3,
                    snr=100,
                    kwargs...)
    maxfcsep = size(cube, 1) ÷ (2 * fwhm) - 1
    # too large separation between companions in the radial patterns
    3 ≤ fc_rad_sep ≤ maxfcsep || error("`fc_rad_sep` should lie ∈[3, $(maxfcsep)], got $fc_rad_sep")

    # compute noise in concentric annuli on the empty frame
    reduced_empty = alg(cube, angles; kwargs...)

    noise, radii = noise_per_annulus(reduced_empty, fwhm, fwhm)
    noise = @view noise[inner_rad:end]
    radii = @view radii[inner_rad:end]
    center_ = center(cube)[2:3]

    psf_size = round(Int, 3 * fwhm)
    iseven(psf_size) && (psf_size += 1)

    psf = construct(psf_model, (psf_size, psf_size))

    angle_branch = 360 / nbranch
    output = zeros(nbranch, length(noise))
    for branch in 0:nbranch, irad in 1:fc_rad_sep
        rads = @view radii[irad:fc_rad_sep:end]
        cube_fake_comps = copy(cube)
        fake_comps = zeros(size(reduced_empty))

        apertures = map(eachindex(rads)) do idx
            A = snr * noise[irad + (idx - 1) * fc_rad_sep]
            injected_flux[idx] = A * psf_flux
            r = rads[idx]
            t = branch * angle_ranch + theta
            inject!(cube_fake_comps, psf, angles; A=A, r=r, theta=t, kwargs...)
            inject!(fake_comps, psf; A=A, r=r, theta=t, kwargs...)
            
            CircularAperture(reverse(r .* sincosd(t) .+ center_), fwhm)
        end
        reduced = alg(cube_fake_comps, angles; kwargs...)
        injected_flux = photometry(apertures, fake_comps).aperture_sum
        recovered_flux = photometry(apertures, reduced .- reduced_empty).aperture_sum

        throughput = recovered_flux ./ injected_flux
        @. throughput[throughput < 0] = 0

        output[branch + 1, irad:fc_rad_sep:end] .= throughput
    end

    return output, radii
end

function noise_per_annulus(frame::AbstractMatrix, separation, fwhm; r0=fwhm)
    cy, cx = center(frame)
    n_annuli = floor(Int, (cy - r0) / separation) - 1
    
    noise_ann = similar(frame, float(eltype(frame)), n_annuli)
    radii = similar(frame, float(eltype(frame)), n_annuli)


    for i in 0:n_annuli-1
        r = r0 + separation * i
        y = cy + r
        
        # find coordinates 
        npoints = floor(Int, 2 * π * r / fwhm)

        apertures = map(range(0, 360, length=npoints)) do theta
            y, x = r .* sincosd(theta)
            y += cy
            x += cx
            CircularAperture(x, y, fwhm/2)
        end
        
        fluxes = photometry(apertures, frame).aperture_sum
        noise_ann[i + 1] = std(fluxes)
        radii[i + 1] = r
    end
    return noise_ann, radii
end

