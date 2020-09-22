# Getting Started

Here is a quick-start guide for people familiar with ADI and experience using tools like [VIP]() or [PyKLIP](). 

## Expected Data Formats

### ADI Cube

For standard ADI data, we store the values in a 3-dimensional array, where the first dimension is temporal, and the remaining dimensions are pixel coordinates. This is how most ADI data are stored on disk (typically in FITS files) and allow specifying operations like a tensor. 

## Algorithms

The following algorithms are implemented:
* [Median Subtraction](@ref med)
* [PCA](@ref pca)
* [NMF](@ref nmf)
* [GreeDS](@ref greeds)

## Processing Patterns

### Full Reduction

Given an algorithm `alg`, we can fully process ADI data by calling `alg` like a function

```julia
julia> alg = PCA(5)

julia> resid = alg(cube, angles)
```

### Reduction Process

The process for producing the flat, residual frame follows this general workflow

1. Create a cube of the speckle approximation, `S`
2. Subtract `S` from the data cube to create the residual cube `R`
3. Derotate `R` frame-by-frame according to the parallactic angle
4. Collapse the derotated `R` 

In ADI.jl this process looks like this:

```julia
using HCIToolbox: collapse, derotate
using ADI: reconstruct, PCA

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
S \approx \mathbf{w} \cdot \mathbf{A}
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
