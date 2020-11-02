#=
# ADI Reduction of $\beta$ Pictoris

This example will walk through a full reduction and some analysis of $\beta$ Pictoris in order to show the basic usage of ADI.jl for high-contrast imaging. $\beta$ Pictoris is a known exoplanet host with many publications regarding the direct images of its substellar companion and circumstellar disk.

What will *not* be covered in this example are the basics of Julia, the fine details of ADI post-processing, or any reference documentation.

---

## Setup

Let's begin by importing the necessary libraries
=#
using ADI
using DataFrames
using HCIDatasets: BetaPictoris
using HCIToolbox
using Plots

## set up plotting
function imshow(img, args...; kwargs...)
    ylim = extrema(axes(img, 1))
    xlim = extrema(axes(img, 2))
    heatmap(img, args...; aspect_ratio=1, xlim=xlim, ylim=ylim, kwargs...)
end;

#=
## Data Reduction

Here we load the data for $\beta$ Pictoris from NaCo at the VLT. You may be prompted to download the data; see [HCIDatasets.jl](https://github.com/JuliaHCI/HCIDatasets.jl) for more details.
=#
cube = BetaPictoris[:cube]
angles = BetaPictoris[:pa];

#=
To reduce the data, we need an *algorithm*. In ADI.jl we currently have median subtraction, PCA, NMF, and fixed-point GreeDS. These algorithms are treated as "objects" in the sense that we initialize them with options and then pass them around inside the ADI.jl framework to retrieve the results we want.

The usage for fitting the speckle estimate, projecting and subtracting this estimate from the target cube, and derotating and collapsing the residual all are encompassed by calling the algorithm as a function.
=#
alg = PCA(10) # 10 components
reduced = alg(cube, angles)
imshow(reduced)

#=
To try out different algorithms, all you should have to do is change one line and re-run the remaining code. Let's briefly explore a few different algorithms
=#
algs = (PCA(10), PCA(20), Median(), GreeDS(10))
reduced_frames = [alg(cube, angles) for alg in algs];
#-
figs = (imshow(reduced, ticks=false) for reduced in reduced_frames)
plot(
    figs...,
    title=["PCA(10)" "PCA(20)" "median subtraction" "GreeDS(10)"],
    layout=(2, 2),
    size=(600, 500)
)

#=
You may want to mask out an interior angle since there is an inner limit for our signal to be a real planet (as opposed to systematics from the optical system or noise). We can mask out an interior circle either before processing with the algorithm or afterwards using `HCIToolbox.mask_circle`.
=#
mask_cube = mask_circle(cube, 10)
mask_reduced = alg(mask_cube, angles)
imshow(mask_reduced)

#=
## S/N and Significance Maps

Now that we have our reduced frame, let's look at the signal-to-noise ratio (SNR, S/N). We use the exact S/N calculation here, implemented in a fast, multi-threaded framework. In order to measure the S/N, though, we need the effective FWHM of our instrument. Normally, we would measure this from an off-axis (or non-coronagraphic) PSF, but for simplicity I'll hard-code a value.
=#
fwhm = 4.6
snrmap = detectionmap(snr, reduced, fwhm)
imshow(snrmap)

#=
If we want to get the statistical significance, we need to convert from the Student-t confidence interval derived in the S/N to the Gaussian significance. We can accomplish this by calling
=#
sigmap = detectionmap(significance, reduced, fwhm)
imshow(sigmap)

#=
Now, lets do some very basic frequentist planet detection by thresholding this significance
=#
sigmap_cutoff = @. sigmap > 5
imshow(sigmap_cutoff)

#=
looks like we've successfully pulled out the companion Beta Pictoris b from the data!

## Contrast Curve

We are also interested in analyzing how the algorithm affects our data, especially calculating the *throughput* and the *contrast curve*. These measure, effectively, how much signal is lost during the subtraction step of the algorithm and give us an idea of what the limits of our algorithm are with our data.

Before we move on, we need to create a PSF model for our data. `HCIToolbox.Kernels` includes some simple functional PSFs or we can use an empirical PSF. We will use the empirical PSF provided by HCIDatasets for our calculations
=#
psf = BetaPictoris[:psf] ./ maximum(BetaPictoris[:psf])
kern_psf = Kernels.Normal(fwhm);
#-
plot(
    imshow(psf),
    imshow(construct(kern_psf, size(psf))),
    layout=2,
    size=(500, 250),
    cbar=false,
    ticks=false
)

#=
and now we can calculate the 5σ contrast curve
=#

cc = contrast_curve(alg, cube, angles, psf; fwhm=fwhm, nbranch=6) |> DataFrame
head(filter(row -> isfinite(row.contrast_corr), cc))

# and lets plot it
plot(
    cc[:distance],
    [cc[:contrast_corr] cc[:contrast]],
    yscale=:log10,
    xlim=(0, NaN),
    label=["Student-t" "Gaussian"],
    ylabel="5-sigma contrast",
    xlabel="radius [px]"
)
#=
You'll notice a pretty severe bump indicating poor contrast where the original companion is! Because this companion is very bright (it can even be seen in the raw frames!) it will bias the contrast measurement.

Typically you'd like to fit the companion signal and remove it in a maximum likelihood framework. For convenience here, let's use the `:cube_empty` entry for `BetaPictoris`, which already has the companion removed.
=#
cube_empty = BetaPictoris[:cube_empty]
reduced_empty = alg(cube_empty, angles)
imshow(reduced_empty)
#-
cc_empty = contrast_curve(alg, cube_empty, angles, psf; fwhm=fwhm, nbranch=6) |> DataFrame
head(filter(row -> isfinite(row.contrast_corr), cc_empty))
#-
plot(
    cc_empty[:distance],
    [cc_empty[:contrast_corr] cc_empty[:contrast]],
    yscale=:log10,
    xlim=(0, NaN),
    label=["Student-t" "Gaussian"],
    ylabel="5-sigma contrast",
    xlabel="radius [px]"
)
