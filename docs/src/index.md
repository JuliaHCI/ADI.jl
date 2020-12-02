```@meta
CurrentModule = ADI
```

# ADI.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/juliahci/ADI.jl)
[![Build Status](https://github.com/juliahci/ADI.jl/workflows/CI/badge.svg?branch=master)](https://github.com/juliahci/ADI.jl/actions)
[![Coverage](https://codecov.io/gh/juliahci/ADI.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliahci/ADI.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

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

## License

This work is distributed under the MIT "expat" license. See [`LICENSE`](https://github.com/juliahci/ADI.jl/blob/master/LICENSE) for more information.
