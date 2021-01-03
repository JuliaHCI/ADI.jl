module ADI

using LinearAlgebra
using ProgressLogging
using Reexport
using Statistics
using StructArrays

@reexport using HCIToolbox

export reconstruct,
       subtract,
       decompose,
       Classic,
       PCA,
       GreeDS,
       TPCA,
       NMF,
       SingleSDI,
       DoubleSDI,
       SliceSDI

"""
    ADI.ADIAlgorithm

This abstract type is used for defining ADI algorithms. Algorithms are stateful objects that define the options for a given algorithm (e.g. the number of components used in PCA decomposition). The most direct usage of an algorithm is to use it to fit HCI data; that is, given a sample of pixels, apply the algorithm using the given options and return an object containing all the necessary information to reconstruct the input data.

See the extended help (`??ADIAlgorithm`) for interface details.

# Extended help
## Interface
To extend `ADIAlgorithm` you may implement the following

| function | default | description |
|----------|---------|-------------|
| `fit(::Alg, data::AbstractMatrix; kwargs...)` | | Fit the data (flattened into a matrix). To support RDI, ensure the `ref` keyword argument is usable (`ref` is also a flattened matrix). |
| `reconstruct(::Alg, cube; kwargs...)` | Uses the reconstruction from the [`ADIDesign`](@ref) fit to the data | Fit the data using the algorithm and return a cube with the estimate of the PSF |
| `subtract(::Alg, cube; kwargs...)` | calls [`reconstruct`](@ref) and subtracts it from `cube` | Fit the data using the algorithm and return a cube that has had the PSF estimate subtracted. |
| `(::ADIAlgorithm)(cube; kwargs...)` | derotates and collapses output from `subtract` | Fully process the data (estimate, subtract, collapse) |

ADI.jl automatically coaxes the `cube` input into a matrix for use with `fit`, appropriately handling the various geometries. If a given algorithm doesn't support the default operations, all that needs to be done is override the default behavior (for an example, see the [`GreeDS`](@ref) implementation).
"""
abstract type ADIAlgorithm end

"""
    ADIDesign


## Interface
To extend `ADIDesign` you may implement the following
| function | default | description |
|----------|---------|-------------|
| `ADI.design(::Design)` | | return the pertinent data to approximate the PSF
| `reconstruct(::Design)` | | return the approximate PSF estimate as a matrix
"""
abstract type ADIDesign end

"""
    LinearDesign <: ADIDesign

A "linear" design implies the use of some linear basis for reconstructing data along with a set of coefficients or weights. Designs which subtype `LinearDesign` must implement the [`ADI.design`](@ref) method, which returns the basis and weights (in that order) that were fit to the data.

`LinearDesign`s have a default implementation of `reconstruct` which uses the matrix product of the weights and the basis, `w * Z`.
"""
abstract type LinearDesign <: ADIDesign end

"""
    ADI.design(::ADIDesign)

Return the pertinent data required to form the PSF approximation. For example, weights and components for PCA/NMF or the median frame for classic ADI.
"""
function design end

"""
    fit(::ADIAlgorithm, cube; [ref], kwargs...)

Given the description of an algorithm and the appropriate options, take the pixels from `cube` and fit them, returning a design ([`ADIDesign`](@ref)) containing the necessary information from the fit (e.g. the principal components from PCA decomposition).

If the algorithm supports reference differential imaging (RDI), the reference cube can be passed by the keyword argument `ref`.
"""
function fit(alg::ADIAlgorithm, cube::AbstractArray{T,3}; kwargs...) where T
    data = flatten(cube)
    if :ref in keys(kwargs)
        kwargs[:ref] isa AbstractArray{<:Any,3} || error("reference data geometry does not match target data")
        ref_data = flatten(kwargs[:ref])
        return fit(alg, data; kwargs..., ref=ref_data)
    end
    return fit(alg, data; kwargs...)
end
function fit(alg::ADIAlgorithm, cube::AnnulusView; kwargs...)
    data = cube()
    if :ref in keys(kwargs)
        kwargs[:ref] isa AnnulusView || error("reference data geometry does not match target data")
        ref_data = kwargs[:ref]()
        return fit(alg, data; kwargs..., ref=ref_data)
    end
    return fit(alg, data; kwargs...)
end

function fit(alg::ADIAlgorithm, cube::MultiAnnulusView; kwargs...)
    if :ref in keys(kwargs)
        kwargs[:ref] isa MultiAnnulusView || error("reference data geometry does not match target data")
        anns = eachannulus(cube)
        ref_anns = eachannulus(kwargs[:ref])
        return StructArray(fit(alg, ann; kwargs..., ref=ref_ann) for (ann, ref_ann) in zip(anns, ref_anns))
    else
        return StructArray(fit(alg, ann; kwargs...) for ann in eachannulus(cube))
    end
end

"""
    reconstruct(::ADIAlgorithm, cube; [ref], kwargs...)

Reconstruct the PSF approximation for the given algorithm, using `ref` as the reference cube if given and supported by the algorithm.

# Examples

```julia
julia> cube, angles = # load data

julia> S = reconstruct(PCA(10), cube);

julia> size(S) == size(cube)
true

julia> flat_res = collapse(cube .- S, angles); # form resid, derotate, and combine
```
"""
function reconstruct(alg::ADIAlgorithm, cube; kwargs...)
    design = fit(alg, cube; kwargs...)
    S = reconstruct(design)
    return expand_geometry(cube, S)
end

# eg MultiAnnulusView
reconstruct(designs::AbstractVector{<:ADIDesign}) = map(reconstruct, designs)

function reconstruct(des::LinearDesign)
    A, weights = design(des)
    return weights * A
end


expand_geometry(::AbstractArray{T,3}, arr) where {T} = expand(arr)
expand_geometry(cube::AnnulusView, arr) = inverse(cube, arr)
expand_geometry(cube::MultiAnnulusView, arrs) = inverse(cube, arrs)

"""
    subtract(::ADIAlgorithm, cube; [ref], kwargs...)

Reconstruct the PSF approximation for the given algorithm and subtract it from `cube`, using `ref` as the reference cube if given and supported by the algorithm.

# Examples

```julia
julia> cube, angles = # load data

julia> R = subtract(PCA(10), cube);

julia> size(R) == size(cube)
true

julia> flat_res = collapse(R, angles); # derotate, and combine
```
"""
function subtract(alg::ADIAlgorithm, cube; kwargs...)
    S = reconstruct(alg, cube; kwargs...)
    return cube .- S
end

"""
    (::ADIAlgorithm)(cube, angles; [ref], kwargs...)

Fully process an ADI data cube using [`reconstruct`](@ref) and collapsing the residuals. Keyword arguments will be passed to `HCIToolbox.collapse!`.
"""
function (alg::ADIAlgorithm)(cube, angles; method=:deweight, kwargs...)
    return collapse!(subtract(alg, cube; kwargs...), angles, method=method)
end

# The core decomposition routines
include("classic.jl")
include("pca.jl")
include("greeds.jl")
include("nmf.jl")
include("sdi.jl")

include("metrics/Metrics.jl")
@reexport using .Metrics

end
