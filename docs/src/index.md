```@meta
CurrentModule = ADI
```

# ADI.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/juliahci/ADI.jl)
[![Build Status](https://github.com/juliahci/ADI.jl/workflows/CI/badge.svg?branch=main)](https://github.com/juliahci/ADI.jl/actions)
[![Coverage](https://codecov.io/gh/juliahci/ADI.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/juliahci/ADI.jl)
[![License](https://img.shields.io/github/license/JuliaHCI/ADI.jl?color=yellow)](https://github.com/JuliaHCI/ADI.jl/blob/main/LICENSE)

[![JOSS](https://joss.theoj.org/papers/32605be405e024fcbd15cd81dfdf9985/status.svg)](https://joss.theoj.org/papers/32605be405e024fcbd15cd81dfdf9985)
[![DOI](https://zenodo.org/badge/250468435.svg)](https://zenodo.org/badge/latestdoi/250468435)

A package for angular differential imaging (ADI) along with its variants, such as reference differential imaging (RDI) and spectral differential imaging (SDI).

## Installation and Setup

ADI.jl is a registered package and can be installed using the Julia package manager. From the Julia REPL, enter Pkg mode (by pressing `]`)

```julia
julia>]

(@v1.5) pkg> add ADI
```

To exit Pkg mode, just backspace. Once the package is installed it can be imported with

```julia
julia> using ADI
```

For more information, see the [Pkg documentation](https://docs.julialang.org/en/v1/stdlib/Pkg/).

## Citations

If you use ADI.jl or derivatives in your work, please consider citing both the JOSS paper and the code record. The JOSS paper citation can be found in [`CITATION.bib`](https://github.com/juliahci/ADI.jl/blob/master/CITATION.bib). The code will have a unique reference for each released version, so visit the [Zenodo record](https://doi.org/10.5281/zenodo.3977789) to grab the BibTeX for whichever version you used.

## Contributing and Support

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

In general contributions should follow [ColPrac](https://github.com/SciML/ColPrac). If you are interested in extending/improving ADI.jl, head to the [discussions](https://github.com/JuliaHCI/ADI.jl/discussions) to reach out. For support with using ADI.jl, please open an [issue](https://github.com/JuliaHCI/ADI.jl/issues/new/) describing the problem and steps to reproduce it.

## License

This work is distributed under the MIT "expat" license. See [`LICENSE`](https://github.com/juliahci/ADI.jl/blob/main/LICENSE) for more information.
