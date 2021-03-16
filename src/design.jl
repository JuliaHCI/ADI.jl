
"""
    ADI.ADIDesign

An abstract type used for holding the output of [`ADI.fit`](@ref). The purpose of these types is to hold the *minimum* information required to reconstruct the PSF estimate (e.g. weights and principal components from PCA). ADI.jl includes two designs- [`ADI.ClassicDesign`](@ref) and [`ADI.LinearDesign`](@ref).

## Interface

When implementing a new algorithm, if your outputs do not fit into either of those designs, you will have to create your own `ADIDesign` with the following methods-

    ADI.design(::Design)

accessor for the pertinent data to approximate the PSF

---

    reconstruct(::Design)

return the approximate PSF estimate as a matrix
"""
abstract type ADIDesign end

# eg MultiAnnulusView
reconstruct(designs::AbstractVector{<:ADIDesign}) = map(reconstruct, designs)

"""
    ADI.LinearDesign(basis, coeffs)

A "linear" design implies the use of some linear basis for reconstructing data along with a set of coefficients or weights. The reconstruction will be the matrix product of the weights and the basis, `w * Z`. 

[`ADI.design`](@ref) will return `(basis, ceoffs)`, and you can also extract them via iteration, like `Z, w = design`.
"""
@concrete struct LinearDesign <: ADIDesign
    basis
    coeffs
end

design(des::LinearDesign) = des.basis, des.coeffs
reconstruct(des::LinearDesign) = des.coeffs * des.basis
Base.iterate(des::LinearDesign, state=1) = iterate(design(des), state)

design(mat::AbstractMatrix) = mat
reconstruct(mat::AbstractMatrix) = mat

"""
    ADI.design(::ADIDesign)

Return the pertinent data required to form the PSF approximation. For example, weights and components for PCA/NMF or the median frame for classic ADI.
"""
function design end
