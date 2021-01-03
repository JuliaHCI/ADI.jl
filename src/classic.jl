
"""
    Classic(;method=median)
    Classic(method)

Classic PSF subtraction using the median of entire data cube. If another statistic is desired, like the `mean`, it can be passed as an argument as long as it supports slicing multi-dimensional arrays.

# References
1. [Marois et al. 2006](http://adsabs.harvard.edu/abs/2006ApJ...641..556M) Angular Differential Imaging: A Powerful High-Contrast Imaging Technique
"""
struct Classic <: ADIAlgorithm
    method
end
Classic(;method=median) = Classic(method)

struct ClassicDesign{FT} <: ADIDesign
    n::Int # number of frames in original data
    frame::FT
end

design(des::ClassicDesign) = des.frame

reconstruct(des::ClassicDesign) = repeat(des.frame, des.n, 1)

function fit(alg::Classic, data::AbstractMatrix; ref=data)
    n = size(data, 1)
    return ClassicDesign(n, alg.method(ref, dims=1))
end
