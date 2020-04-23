module ADI

using HCIToolbox
import MultivariateStats
using Statistics

const mvs = MultivariateStats

export pca, transform, reconstruct

abstract type ADIDesign end

"""
    transform(::ADIDesign, cube)
    transform(::ADIDesign, matrix)
"""
transform(d::ADIDesign, cube::AbstractArray{T, 3}) where T = transform(d, flatten(cube))

"""
    reconstruct(::ADIDesign, weights)
"""
reconstruct


"""
    reduce(::ADIDesign, [cube, angles]; method=median, deweight=true, fill=0)

Reduces an ADI Design matrix by computing the residual cube and collapsing it. The keyword arguments will be passed to `HCIToolbox.collapse!`. If `cube` and `angles` are provided, they will be used for calculating the residual. Otherwise, the cube and angles used to create the design will be used.
"""
function Base.reduce(d::ADIDesign, cube::AbstractArray{T,3}, angles::AbstractVector; kwargs...) where T
    w = transform(d, cube)
    S = reshape(reconstruct(d, w), size(cube))
    R = _cube .- S
    return collapse!(R, angles; kwargs...)
end

Base.reduce(d::ADIDesign; kwargs...) = collapse!(d._cube .- d.S, d._angs; kwargs...)


# The core decomposition routines
include("pca.jl")

end
