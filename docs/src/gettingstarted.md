# [Getting Started](@id gettingstarted)

Here is a quick-start guide for people familiar with ADI and experience using tools like [VIP](https://github.com/vortex-exoplanet/VIP) or [PyKLIP](https://pyklip.readthedocs.io/en/latest/). For installation and setup information, see the [Installation and Setup](@ref) section.

## Expected Data Formats

### ADI Cube

For standard ADI data, we store the values in a 3-dimensional array, where the first dimension is temporal, and the remaining dimensions are pixel coordinates. This is how most ADI data are stored on disk (typically in FITS files) and allow specifying operations like a tensor.

### SDI Cube/Tensor

For standard SDI data, we store the values in a 4-dimensional array, where the first dimension is spectral, the second is temporal, and the remaining dimensions are pixel coordinates. This is how most SDI data are stored on disk (typically in FITS files) and allow specifying operations like a tensor. In addition to the SDI tensor and parallactic angles, the list of wavelengths are required (for scaling speckles) and a spectral template can be used.

## Algorithms

The following algorithms are implemented:
* [Median Subtraction](@ref med)
* [PCA](@ref pca)
* [NMF](@ref nmf)
* [GreeDS](@ref greeds)

## Processing Patterns

### Full Frame ADI Reduction

Given an algorithm `alg`, we can fully process ADI data by calling `alg` like a function

```julia
julia> alg = PCA(ncomps=5)

julia> resid = alg(cube, angles)
```

### Full Frame RDI Reduction

The only difference here is the inclusion of a reference cube.

```julia
julia> alg = PCA(ncomps=5)

julia> resid = alg(cube, angles, cube_ref)
```

### Reduction Process

The process for producing the flat, residual frame follows this general workflow

1. Create a cube of the speckle approximation, `S`
2. Subtract `S` from the data cube to create the residual cube `R`
3. Derotate `R` frame-by-frame according to the parallactic angle
4. Collapse the derotated `R`

In ADI.jl this process looks like this:

```julia
using ADI

cube, angles = # load data
S = reconstruct(PCA(10), cube, angles)
R = cube .- S
R_derotate = derotate(R, angles)
resid = collapse(R_derotate)
# or, more succinctly
resid = collapse(R, angles)
```

Notice how the only part of this specific to the algorithm is [`reconstruct`](@ref)? This lets us have the super-compact functional form from above without having to copy the common code for each algorithm.

### Decomposition

For certain types of ADI algorithms, a convenient linear form is used for the speckle approximation
```math
\mathbf{S} \approx \mathbf{w} \cdot \mathbf{A}
```
Algorithms which share this attribute share the abstract type `ADI.LinearAlgorithm`, and we can retrieve these two matrices via [`decompose`](@ref).

```julia
using ADI: decompose, reconstruct, PCA
cube, angles = # load data
A, w = decompose(PCA(10), cube, angles)
S = reconstruct(PCA(10), A, w)
S == w * A
# output
true
```

## Comparison to VIP

ADI.jl took a lot of ideas from VIP and expanded them using the power of Julia. To begin with, Julia typically has smaller and more self-contained packages, so most of the basic image-processing that is used here is actually written in the [HCIToolbox.jl](https://github.com/JuliaHCI/HCIToolbox.jl) package. In the future, I have plans to incorporate forward-modeling distributions in [Firefly.jl](https://github.com/JuliaHCI/Firefly.jl), which currently is an active research project.

Some technical distinctions to VIP
* Julia is 1-indexed. This means all the positions for apertures, bounds, images, etc. start at 1. This is distinct from 0-based indexing in python, but is equivalent to the indexing in DS9 and IRAF.
* Julia's `std` uses the sample statistic ($n-1$ degrees of freedom) while numpy's `std` uses the population statistic ($n$ degrees of freedom). This may cause very slight differences in measurements that rely on this.
* Aperture mapping - many of the [`Metrics`](@ref) are derived by measuring statistics in an annulus of apertures. In VIP, this ring is not equally distributed- the angle between apertures is based on the exact number of apertures rather than the integral number of apertures that are actually measured. In ADI.jl the angle between apertures is evenly distributed. The same number of pixels are discarded in both packages, but in VIP they all end up in the same region of the image (see [this figure](assets/aperture_masks.png)).
* Collapsing - by default VIP collapses a cube by derotating it then finding the median frame. In ADI.jl, the default collapse method is a weighted sum using the inverse of the temporal variance for weighting. This is documented in `HCIToolbox.collapse` and can be overridden by passing the keyword argument `method=median` or whichever statistical funciton you want to use.

The biggest difference, though, is Julia's multiple-dispatch system and how that allows ADI.jl to *do more with less code*. For example, the [`GreeDS`](@ref) algorithm was designed explicitly for [`PCA`](@ref), but the formalism of it is more generic than that. Rather than hard-coding in PCA, the GreeDS algorithm was written generically, and Julia's multiple-dispatch  allows the use of, say, [`NMF`](@ref) instead of PCA. By making the code *generic* and *modular*, ADI.jl enables rapid experimentation with different post-processing algorithms and techniques as well as minimizing the code required to implement a new algorithm and be able to fully use the ADI.jl API.

