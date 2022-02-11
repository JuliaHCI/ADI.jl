# ADI.jl

[![Build Status](https://github.com/juliahci/ADI.jl/workflows/CI/badge.svg?branch=main)](https://github.com/juliahci/ADI.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/A/ADI.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/juliahci/ADI.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/juliahci/ADI.jl)
[![License](https://img.shields.io/github/license/JuliaHCI/ADI.jl?color=yellow)](LICENSE)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliahci.github.io/ADI.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliahci.github.io/ADI.jl/dev)
[![JOSS](https://joss.theoj.org/papers/32605be405e024fcbd15cd81dfdf9985/status.svg)](https://joss.theoj.org/papers/32605be405e024fcbd15cd81dfdf9985)
[![DOI](https://zenodo.org/badge/250468435.svg)](https://zenodo.org/badge/latestdoi/250468435)

A package for angular differential imaging (ADI) post-processing algorithms.

## Installation

ADI.jl is a registered package and can be installed using the Julia package manager. From the Julia REPL, enter Pkg mode (by pressing `]`)

```julia
julia>]

(@v1.5) pkg> add ADI
```

To exit Pkg mode, just backspace. Once the package is installed it can be imported with

```julia
julia> using ADI
```

To exit Pkg mode, just backspace. Once the package is installed it can be imported with
For more information, see the [Pkg documentation](https://docs.julialang.org/en/v1/stdlib/Pkg/).

## Citations

If you use ADI.jl or derivatives in your work, please consider citing both the JOSS paper and the code record. The JOSS paper citation can be found in [`CITATION.bib`](CITATION.bib). The code will have a unique reference for each released version, so visit the [Zenodo record](https://doi.org/10.5281/zenodo.3977789) to grab the BibTeX for whichever version you used.

## Usage

The following is an extremely brief PCA reduction of an ADI cube. Please see the [documentation](https://juliahci.github.io/ADI.jl/dev/) for further usage, tutorials, and api reference.

```julia
julia> using ADI

julia> cube, angles = # load data

julia> alg = PCA(ncomps=10)

julia> flat_resid = alg(cube, angles) # ADI

julia> flat_resid_rdi = alg(cube, angles; ref=cube_ref) # flexible RDI
```

get the S/N and significance

```julia
julia> fwhm = # PSF fwhm in pixels

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

julia> first(df, 5)
```

## Contributing and Support

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

In general contributions should follow [ColPrac](https://github.com/SciML/ColPrac). If you are interested in extending/improving ADI.jl, head to the [discussions](https://github.com/JuliaHCI/ADI.jl/discussions) to reach out. For support with using ADI.jl, please open an [issue](https://github.com/JuliaHCI/ADI.jl/issues/new/) describing the problem and steps to reproduce it.

## License

This package is licensed under the MIT Expat license. See [LICENSE](LICENSE) for more information.

---

**Author's Note**: This package is still under active development and is subject to change. Anything from minor behind-the-scenes details to large-scale design can change as I incorporate more methods into ADI.jl. I don't plan on spending much time with deprecation warnings throughout this process, since that limits my ability to experiment with implementation ideas and design goals. This package follows [semantic versioning](https://semver.org/), so an upgrade from `0.6` to `0.7` may be breaking and I recommend anybody using this package to browse the release notes for changes. Once ADI.jl is somewhat stable, I'll release a version `1.0`, at which point I'll worry about deprecations and other long-term usability considerations.
