module ADI

using HCIToolbox
using Statistics

export pca, pairet, transform, reconstruct, projection

abstract type ADIDesign{T<:AbstractArray, V<:AbstractVector} end

"""
    reconstruct(::ADIDesign, cube) -> cube
    reconstruct(::ADIDesign, matrix) -> matrix

Reconstrucst an approximation of the input using the design.

# See Also
[`reduce`](@ref)
"""
reconstruct(d::ADIDesign, cube::AbstractArray{T, 3}) where T = reconstruct(d, flatten(cube)) |> expand

"""
    reduce(::ADIDesign; method=median, deweight=true, fill=0) -> matrix
    reduce(::ADIDesign, cube, [angles]; method=median, deweight=true, fill=0) -> matrix

Reduces an ADI Design matrix by computing the residual of the reconstructed cube and the target cube, then collapsing it. The keyword arguments will be passed to `HCIToolbox.collapse!`. 

If `cube` and `angles` are provided, referential differential imaging (RDI) will be done by reconstructing `cube` using the input design. If `angles` are not provided, the same angles used for the design construction will be used.

# See Also
[`reconstruct`](@ref)
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
