module ADI

using HCIToolbox
using Statistics
using Parameters
using UnPack

export pca, pairet, transform, reconstruct, projection

export ADIAlgorithm,
       design,
       Median

"""
    ADIAlgorithm <: Function
    
This abstract type is used for defining ADI algorithms. See the extended help (`??ADIAlgorithm`) for interface details

# Extended Help
## Interface
To extend `ADIAlgorithm` you may implement the following
| function | default | description |
|----------|---------|-------------|
| `ADI.design`
"""
abstract type ADIAlgorithm <: Function end

"""
    design(::ADIAlgorithm)(cube, angles, args...; kwargs...)

Return some structure corresponding to the given algorithm. By default, will return a `NamedTuple` of the `cube`, `angles`, and the output of [`(::ADIAlgorithm)`](@ref)
"""
function design(alg::ADIAlgorithm, cube, angles, args...; kwargs...)
    output = alg(cube, angles, args...; kwargs...)
    return (cube=cube, angles=angles, reduced=output)
end

include("median.jl")

# # The core decomposition routines
# include("pca.jl")
# include("pairet.jl")

# include("medsub.jl")


# abstract type ADIDesign{T,A<:AbstractArray{T},V<:AbstractVector{T}} end

# """
#     reconstruct(::ADIDesign, cube) -> cube
#     reconstruct(::ADIDesign, matrix) -> matrix

# Reconstrucst an approximation of the input using the design.

# # See Also
# [`reduce`](@ref)
# """
# reconstruct(d::ADIDesign, cube::AbstractArray{T, 3}) where T = 
#     reconstruct(d, flatten(cube)) |> expand

# """
#     reduce(::ADIDesign; method=:deweight, fill=0, degree=Linear()) -> matrix
#     reduce(::ADIDesign, cube, [angles]; method=:deweight, fill=0, degree=Linear()) -> matrix

# Reduces an ADI Design matrix by computing the residual of the reconstructed cube and the target cube, then collapsing it. The keyword arguments will be passed to `HCIToolbox.collapse!`. 

# If `cube` and `angles` are provided, referential differential imaging (RDI) will be done by reconstructing `cube` using the input design. If `angles` are not provided, the same angles used for the design construction will be used.

# # See Also
# [`reconstruct`](@ref)
# """
# function Base.reduce(d::ADIDesign, cube::AbstractArray{T,3}, angles::AbstractVector=d.angles; kwargs...) where T
#     reconstructed = reconstruct(d, cube)
#     R = cube .- reconstructed
#     return collapse!(R, angles; kwargs...)
# end

# Base.reduce(d::ADIDesign; kwargs...) = collapse(d.S, d.angles; kwargs...)

using Reexport

include("metrics/Metrics.jl")
@reexport using .Metrics

end
