using Statistics

"""
    contrast_curve
"""
function contrast_curve(design::ALG,
                        cube=design.cube,
                        angles=design.angles,
                        psf;
                        fwhm,
                        pxscale,
                        starphot=1,
                        sigma=5,
                        nbranch=1,
                        theta=0,
                        inner_rad=1,
                        fc_snr=100,
                        transmission=nothing) where {ALG<:ADIDesign}



end

"""
    throughput
"""
function throughput(alg,
                    cube::AbstractArray{T,3}=design.cube,
                    angles=design.angles,
                    psf;
                    fwhm,
                    pxscale,
                    sigma=5,
                    nbranch=1,
                    theta=0,
                    inner_rad::Integer=1,
                    fc_rad_sep=3,
                    fc_snr=100,
                    kwargs...) where {T,ALG<:ADIDesign}
    maxfcsep = size(cube, 1) ÷ (2 * fwhm) - 1
    # too large separation between companions in the radial patterns
    3 ≤ fc_rad_sep ≤ maxfcsep || error("`fc_rad_sep` should lie ∈[3, $(maxfcsep)], got $fc_rad_sep")

    Γ = median(fwhm)

    # compute noise in concentric annuli on the empty frame
    frame_empty = ALG()
    
end

