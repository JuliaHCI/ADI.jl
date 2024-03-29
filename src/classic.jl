
"""
    Classic(;method=median)
    Classic(method)

Classic PSF subtraction using the median of entire data cube. If another statistic is desired, like the `mean`, it can be passed as an argument as long as it supports slicing multi-dimensional arrays.

# References
1. [Marois et al. 2006](http://adsabs.harvard.edu/abs/2006ApJ...641..556M) Angular Differential Imaging: A Powerful High-Contrast Imaging Technique
"""
@concrete struct Classic <: ADIAlgorithm
    method
end
Classic(;method=median) = Classic(method)

"""
    ADI.ClassicDesign(n, frame)

Output for the [`Classic`](@ref) algorithm which contains the static frame unrolled into a vector (with size `Npx`). [`reconstruct`](@ref) will tile this vector in a non-allocating way `n` times to appear like a flattened cube. [`ADI.design`](@ref) will return the static frame .
"""
struct ClassicDesign{FT} <: ADIDesign
    n::Int # number of frames in original data
    frame::FT
end

design(des::ClassicDesign) = des.frame

reconstruct(des::ClassicDesign) = TiledArray(des.frame, des.n)

function fit(alg::Classic, data::AbstractMatrix; ref=data, kwargs...)
    n = size(data, 2)
    return ClassicDesign(n, alg.method(ref, dims=2))
end

#######################################

struct TiledArray{T,N,AT<:AbstractArray{T,N},RT} <: AbstractArray{T,N}
    parent::AT
    tilerange::RT
end

Base.size(t::TiledArray) = map(length, axes(t))
Base.axes(t::TiledArray) = (Base.front(axes(t.parent))..., t.tilerange)
Base.parent(t::TiledArray) = t.parent

TiledArray(parent::AbstractArray, n::Int) = TiledArray(parent, Base.OneTo(n))

Base.@propagate_inbounds function Base.getindex(t::TiledArray{T,N}, idxs::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(t, idxs...)
    return t.parent[Base.front(idxs)..., firstindex(t.parent, 3)]
end
