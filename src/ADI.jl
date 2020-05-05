module ADI

using HCIToolbox
using Statistics

export pca, pairet, transform, reconstruct, projection

abstract type ADIDesign{T<:AbstractArray, V<:AbstractVector} end

"""
    reconstruct(::ADIDesign, cube) -> cube
    reconstruct(::ADIDesign, matrix) -> matrix
"""
reconstruct(d::ADIDesign, cube::AbstractArray{T, 3}) where T = reconstruct(d, flatten(cube)) |> expand

"""
    reduce(::ADIDesign; method=median, deweight=true, fill=0)
    reduce(::ADIDesign, cube, [angles]; method=median, deweight=true, fill=0)

Reduces an ADI Design matrix by computing the residual cube and collapsing it. The keyword arguments will be passed to `HCIToolbox.collapse!`. If `cube` and `angles` are provided, they will be used for calculating the residual. Otherwise, the cube and angles used to create the design will be used.
"""
function Base.reduce(d::ADIDesign, cube::AbstractArray{T,3}, angles::AbstractVector=d.angles; kwargs...) where T
    reconstructed = reconstruct(d, cube)
    R = cube .- reconstructed
    return collapse!(R, angles; kwargs...)
end

Base.reduce(d::ADIDesign; kwargs...) = collapse(d.S, d.angles; kwargs...)

# The core decomposition routines
include("pca.jl")
include("pairet.jl")

end
