#=
# [Exploring spatial and temporal filtering](@id ex2)

This example will walk through a more advanced reduction exploiting geometric and temporal filtering techniques.

---

## Setup

Let's begin by importing the necessary libraries. You may need to add these packages if they are not already on your system.
```julia
(@v1.5) pkg> add HCIDatasets Plots
```
In addition, a `Project.toml` file exists in the [examples/](https://github.com/JuliaHCI/ADI.jl/tree/main/examples) folder with the necessary dependencies. Start the REPL in the base ADI.jl folder, then from Pkg mode
```julia
(@v1.5) pkg> activate examples
(examples) pkg> instantiate
julia> include("examples/geometries.jl")
```
---
=#
using ADI
using HCIDatasets: HR8799
using Plots

## set up plotting
function imshow(img; kwargs...)
    xs, ys = axes(img)
    heatmap(xs, ys, transpose(img); aspect_ratio=1, 
            xlim=extrema(xs), ylim=extrema(ys), kwargs...)
end;
#=
## Initial Reduction

First, lets load the data and do some initial reductions
=#
cube, angles = HR8799[:cube, :pa]
fwhm = 8.2

name = "PCA(10) - Full Frame"
res = PCA(10)(cube, angles)
results = Dict(name => res)
imshow(res, title=name)
#=
The four companions (HR8799b,c,d, and e) should be evident in this residual.

## Geometric Filtering

Geometric filtering is the process of using sub-region(s) in each frame to use as the target and reference libraries. Geometric filtering is useful to select pixels with similar noise distributions, such as an annulus, as well as to reduce the total number of pixels, increasing runtime performance.

For our example, lets start by looking at annuli for each companion and optimizing the number of principal components. We can use [`AnnulusView`](@ref) to wrap `cube` and filter the input. [`AnnulusView`](@ref) works by calculating the *indices* for a single annulus and creates a view into `cube` with those indices. It does not copy any data, and lets us access only the pixels within the annulus.
=#
av = AnnulusView(cube; inner=31, outer=51)
imshow(av[:, :, 1], xlim=(198, 304), ylim=(198, 304))
#=
because this annulus should have similar noise characteristics, we can automatically choose the number of components by trying to minimize the noise in the residual frame. The algorithm will increase `ncomps` until the noise does not improve by at least `noise_error`.
=#
alg = PCA(:noise, noise_error=25)
res_e = alg(av, angles)
imshow(res_e, title="PCA - HR8799e", xlim=(198, 304), ylim=(198, 304))
#
av = AnnulusView(cube; inner=57, outer=77)
res_d = alg(av, angles)
imshow(res_d, title="PCA - HR8799d", xlim=(170, 332), ylim=(170, 332))
#
av = AnnulusView(cube; inner=86, outer=106)
res_c = alg(av, angles)
imshow(res_c, title="PCA - HR8799c", xlim=(142, 360), ylim=(142, 360))
#
av = AnnulusView(cube; inner=163, outer=183)
res_b = alg(av, angles)
imshow(res_b, title="PCA - HR8799b", xlim=(57, 445), ylim=(57, 445))
#
combined_res = sum([res_e, res_d, res_c, res_b])
imshow(combined_res, title="PCA - Annular", xlim=(57, 445), ylim=(57, 445))
#=
## Handling multiple annuli

Instead of doing each annulus by hand, we can use [`MultiAnnulusView`](@ref) to assemble the annuli
=#
width = 20
radii = [41, 67, 96, 173]
mav = MultiAnnulusView(cube, width, radii);
#
res = alg(mav, angles)

name = "PCA(noise_error=25) - Annular"
push!(results, name => res)
imshow(res, title=name, xlim=(57, 445), ylim=(57, 445))
#=
When reducing [`MultiAnnulusView`](@ref), we can specify a separate algorithm for each annulus, and even combine that with [`Framewise`](@ref)
=#
algs = [PCA(6), PCA(5), PCA(5), PCA(3)]
res = process(algs, mav, angles)

name = "PCA([6, 5, 5, 3]) - Annular"
push!(results, name => res)
imshow(res, title=name, xlim=(57, 445), ylim=(57, 445))
#=
## Framewise reduction and temporal filtering

As pointed out in [Marois et al. 2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...641..556M/abstract), we can filter the reference data in the time domain by rejecting frames which have not rotated a sufficient amount. In order to enable this, we have to fit each frame in the target data cube separately, since we calculate the frames to reject for each of them. We call this process *framewise* reduction, which is enabled by wrapping an algorithm in [`Framewise`](@ref).

By default, [`Framewise`](@ref) will reject frames which have not rotated at least 1 FWHM, set by the `delta_rot` keyword argument.
=#
fr_alg = Framewise(algs, delta_rot=1)
res = fr_alg(mav, angles; fwhm=fwhm)

name = "PCA([6, 5, 5, 3]) - Annular\n(delta_rot=1)"
push!(results, name => res)
imshow(res, title=name, xlim=(57, 445), ylim=(57, 445))
#=
According to [Absil et al. 2013](https://ui.adsabs.harvard.edu/abs/2013A%26A...559L..12A/abstract), a slightly better contrast can be reached for the innermost annuli if we consider a `delta_rot` as small as 0.1 FWHM. This is because at very small separation, the effect of speckle correlation is more significant than self-subtraction.

We can set `delta_rot` as a tuple or vector to use varying parallactic angle thresholds for each annulus.
=#
fr_alg = Framewise(algs, delta_rot=(0.1, 1))
res = fr_alg(mav, angles; fwhm=fwhm)

name = "PCA([6, 5, 5, 3]) - Annular\n(delta_rot=(0.1, 1))"
push!(results, name => res)
imshow(res, title=name, xlim=(57, 445), ylim=(57, 445))
#=
We can also limit the number of frames used in the reference library, I.E. the `k` nearest frames. For example, to recreate the algorithm derived in § 5.2 of [Marois et al. 2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...641..556M/abstract), which uses the 4 closest frames (which have rotated at least 1.5 FWHM) to construct the reference library.
=#
## subtract median frame
resid_cube = subtract(Classic(), cube)
## form annuli
mav = MultiAnnulusView(resid_cube, width, radii)
## 1.5 FWHM rotation, keep closest 4 frames
fr_alg = Framewise(Classic(), delta_rot=1.5, limit=4)
res = fr_alg(mav, angles; fwhm=fwhm)

name = "Classic - Annular\n(delta_rot=1.5, limit=4)"
push!(results, name => res)
imshow(res, title=name, xlim=(57, 445), ylim=(57, 445))
#=
## Gallery of results
=#
p = plot(layout=(2, 3), ticks=false, xlim=(57, 445), ylim=(57, 445),
        aspect_ratio=1, size=(800, 400), titlefontsize=9)
for (i, (name, res)) in enumerate(results)
    heatmap!(transpose(res); title=name, sp=i)
end
p #hide
#=
### S/N maps
For comparing reductions, we'll use signal-to-noise ratio (S/N, SNR) of the residual frame. All the frames below are shown on the same color scale, from S/N=0 to 13.
=#
p = plot(layout=(2, 3), size=(800, 550), clim=(0, 13),
        cbar=false, ticks=false, aspect_ratio=1, xlim=(57, 445),
        ylim=(57, 445), titlefontsize=9)
for (i, (name, res)) in enumerate(results)
    sn = detectionmap(res, fwhm)
    heatmap!(transpose(sn); title=name, sp=i)
end
p #hide