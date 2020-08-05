module ADI

using HCIToolbox
using Statistics
using Parameters
using UnPack

# export pca, pairet, transform, reconstruct, projection

export reconstruct,
       decompose,
       Median,
       PCA,
       Pairet

"""
    ADI.ADIAlgorithm <: Function

This abstract type is used for defining ADI algorithms. See the extended help (`??ADIAlgorithm`) for interface details.

# Extended help
## Interface
To extend `ADIAlgorithm` you may implement the following
| function | default | description |
|----------|---------|-------------|
| `ADI.reconstruct` | | Subroutine for creating the full reconstructed cube with the PSF |
| `(::ADIAlgorithm)` | subtracts output of `reconstruct`, then derotates and collapses | Subroutine for returning the reduced residual cube |
"""
abstract type ADIAlgorithm <: Function end

"""
    reconstruct(::ADIAlgorithm, cube, angles, [cube_ref]; kwargs...)

Reconstruct the PSF approximation for the given algorithm, using `cube_ref` as the reference cube if given.
"""
function reconstruct end

"""
    (::ADIAlgorithm)(cube, angles, [cube_ref]; kwargs...)

Fully process an ADI data cube using [`reconstruct`](@ref) and collapsing the residuals. Keyword arguments will be passed to `HCIToolbox.collapse!`.
"""
function (alg::ADIAlgorithm)(cube, angles, args...; kwargs...)
    S = reconstruct(alg, cube, angles, args...; kwargs...)
    residual_cube = cube .- S
    return collapse!(residual_cube, angles; kwargs...)
end

include("median.jl")

"""
    ADI.LinearAlgorithm <: ADI.ADIAlgorithm

This abstract type is used for defining linear ADI algorithms. See the extended help (`??LinearAlgorithm`) for interface details.

# Extended help
## Interface
To extend `LinearAlgorithm` you may implement the following

| function | default | description |
|----------|---------|-------------|
| `ADI.fit` | | Subroutine for fitting the linear basis and coefficients as unrolled matrices |
| `ADI.reconstruct` | Computes the inner product of the design matrix and weights from `decompose` | Subroutine for creating the full reconstructed cube with the PSF |
| `(::LinearAlgorithm)` | subtracts output of `reconstruct`, then derotates and collapses | Subroutine for returning the reduced residual cube |
"""
abstract type LinearAlgorithm <: ADIAlgorithm end

"""
    ADI.decompose(::LinearAlgorithm, cube, angles, [cube_ref]; kwargs...)
"""
function decompose end

function reconstruct(alg::LinearAlgorithm, cube, angles, args...; kwargs...)
    # assumed sizes are (n, Npx) (n, M)
    basis, weights = decompose(alg, cube, angles, args...; kwargs...)
    return weights' * basis |> expand
end

# The core decomposition routines
include("pca.jl")
include("pairet.jl")

using Reexport

include("metrics/Metrics.jl")
@reexport using .Metrics

end
