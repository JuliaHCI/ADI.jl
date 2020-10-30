# ADI.jl

[![Build Status](https://github.com/juliahci/ADI.jl/workflows/CI/badge.svg?branch=master)](https://github.com/juliahci/ADI.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/A/ADI.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/juliahci/ADI.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliahci/ADI.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliahci.github.io/ADI.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliahci.github.io/ADI.jl/dev)
[![DOI](https://zenodo.org/badge/250468435.svg)](https://zenodo.org/badge/latestdoi/250468435)

A package for angular differential imaging (ADI) post-processing algorithms.

## Installation

From the Julia REPL

```julia
julia> ]

(@v1.5) pkg> add ADI
```

## Usage

The following is an extremely brief PCA reduction of an ADI cube. Please see the [documentation](https://juliahci.github.io/ADI.jl/dev/) for further usage, tutorials, and api reference.

```julia
julia> using ADI

julia> cube, angles = # load data

julia> alg = PCA(ncomps=10)

julia> flat_resid = alg(cube, angles) # ADI

julia> flat_resid_rdi = alg(cube, angles, cube_ref) # flexible RDI
```

get the S/N and significance

```julia
julia> fwhm = # PSF fswhm in pixels

julia> snmap = detectionmap(snr, flat_residual, fwhm)

julia> sigmap = detectionmap(significance, flat_residual, fwhm)
```

get the contrast curve

```julia
julia> psf = # load psf or choose from HCIToolbox.Kernels

julia> cc = contrast_curve(alg, cube, angles, psf; fwhm=fwhm)
```

which can be easily loaded into a `DataFrame` or anything adopting the Tables.jl interface.

```julia
julia> using DataFrames

julia> df = DataFrame(cc)

julia> head(df)
```

## License

This package is licensed under the MIT Expat license. See [LICENSE](LICENSE) for more information.
