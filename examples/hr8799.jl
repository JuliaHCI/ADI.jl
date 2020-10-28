#=
# ADI Reduction of HR8799

This example will walk through a full reduction and some analysis of HR8799 in order to show the basic usage of ADI.jl for high-contrast imaging. HR8799 is a known exoplanet host with many publications regarding the direct images of 4 sub-stellar companions.

What will *not* be covered in this example are the basics of Julia, the fine details of ADI post-processing, or any reference documentation.

---

## Setup

Let's begin by importing the necessary libraries
=#
using ADI
using DataFrames
using HCIDatasets: HR8799
using HCIToolbox
using Plots

## set up plotting
imshow(args...; kwargs...) = heatmap(args...; aspect_ratio=1, xlims=(1, 501), ylims=(1, 501), kwargs...)

#=
## Data Reduction

Here we load the data for HR8799 from Keck NIRC2/Vortex Coronagraph. You may be prompted to download the data; see HCIDatasets.jl for more details.
=#
cube = HR8799[:cube]
angles = HR8799[:pa];

#=
To reduce the data, we need an *algorithm*. In ADI.jl we currently have median subtraction, PCA, NMF, and fixed-point GreeDS. These algorithms are treated as "objects" in the sense that we initialize them with options and then pass them around inside the ADI.jl framework to retrieve the results we want.

The usage for fitting the speckle estimate, projecting and subtracting this estimate from the target cube, and derotating and collapsing the residual all are encompassed by calling the algorithm as a function. To try out different algorithms, all you should have to do is change this one line and re-run the remaining code.
=#
alg = PCA(10) # 10 components
reduced = alg(cube, angles)
imshow(reduced)

#=
You may want to mask out an interior angle since there is an inner limit for our signal to be a real planet (as opposed to systematics from the optical system or noise). We can mask out an interior circle either before processing with the algorithm or afterwards using `HCIToolbox.mask_circle`.
=#
mask_cube = mask_circle(cube, 16)
mask_reduced = alg(mask_cube, angles)
imshow(mask_reduced)

#=
## S/N and Significance Maps

Now that we have our reduced frame, let's look at the signal-to-noise ratio (SNR, S/N). We use the exact S/N calculation here, implemented in a fast, multi-threaded framework. In order to measure the S/N, though, we need the effective FWHM of our instrument. Normally, we would measure this from an off-axis (or non-coronagraphic) PSF, but for simplicity I'll hard-code a value.
=#
fwhm = 8.637
snrmap = detectionmap(snr, reduced, fwhm)
imshow(snrmap)

#=
If we want to get the statistical significance, we need to convert from the student-t confidence interval derived in the S/N to the Gaussian significance. We can accomplish this by calling
=#
sigmap = detectionmap(significance, reduced, fwhm)
imshow(sigmap)

#=
Now, lets do some very basic frequentist planet detection by thresholding this significance
=#
sigmap_cutoff = @. ifelse(sigmap > 5, sigmap, NaN)
imshow(sigmap_cutoff)

#=
looks like we've successfully pulled out HR8799b,c,d, and e from the data!

## Contrast Curve

We are also interested in analyzing how the algorithm affects our data, especially calculating the *throughput* and the *contrast curve*. These measure, effectively, how much signal is lost during the subtraction step of the algorithm and give us an idea of what the limits of our algorithm are with our data.

Before we move on, we need to create a PSF model for our data. `HCIToolbox.Kernels` includes some simple functional PSFs in the absence of an empirical PSF.
=#
psf = Kernels.Normal(fwhm);

#=
and now we can calculate the 5σ contrast curve
=#

cc = contrast_curve(alg, cube, angles, psf; fwhm=fwhm, nbranch=3) |> DataFrame
head(cc)

# and lets plot it to see how we perform
plot(
    cc[:distance],
    [cc[:contrast_corr] cc[:contrast]],
    ls=[:solid :dash],
    c=1,
    label=["Student-t" "Gaussian"],
    ylabel="5-sigma contrast",
    xlabel="radius [px]"
)
#=
Hmm, there's some pecularities! You'll notice pretty sever bumps indicating poor contrast in annuli where the 4 companions are! Because these companions are pretty significant (statistically) they will bias the contrast measurement.

Typically you'd like to fit the companion signal and remove it in a maximum likelihood framework. For convenience here, though, I am going to create a cube from the speckle estimate, which should be free from companion signal. This is not a rigorous alternative, though, since this cube's noise will be whitened by the PCA process. It will be good enough to demonstrate the code, though.
=#
no_comp_cube = reconstruct(alg, cube, angles)
imshow(alg(no_comp_cube, angles))
#-
cc_no_comp = contrast_curve(alg, no_comp_cube, angles, psf; fwhm=fwhm, nbranch=3) |> DataFrame
head(cc_no_comp)
#-
plot(
    cc_no_comp[:distance],
    [cc_no_comp[:contrast_corr] cc_no_comp[:contrast]],
    ls=[:solid :dash],
    c=1,
    label=["Student-t" "Gaussian"],
    ylabel="5-sigma contrast",
    xlabel="radius [px]"
)
