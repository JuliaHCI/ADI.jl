using Statistics

struct PCADesign{T<:AbstractArray, V<:AbstractVector} <: ADIDesign{T, V}
    pca::mvs.PCA
    S::T
    _cube::T
    _angs::V
end

function Base.show(io::IO, d::PCADesign{T}) where T
    n = mvs.indim(d.pca)
    p = mvs.outdim(d.pca)
    pr = mvs.principalratio(d.pca)
    print("PCADesign{$T}(ncomps=$p, D=$n, pratio=$pr)")
end

transform(d::PCADesign, arr::AbstractMatrix) = mvs.transform(d.pca, arr)
reconstruct(d::PCADesign, w::AbstractMatrix) = mvs.reconstruct(d.pca, w)


"""
    pca(cube, angles, ncomps; pratio=0.99)

Use principal component analysis (PCA) to reduce data cube.
"""
function pca(cube::AbstractArray, angles::AbstractVector; ncomps, pratio = 1, mean=0)
    # transform cube
    nf, ny, nx = size(cube)
    X = reshape(cube, nf, ny * nx)
    P = mvs.fit(mvs.PCA, X; maxoutdim = ncomps, pratio = pratio, mean = mean)
    w = mvs.transform(P, X)
    S = mvs.reconstruct(P, w)
    size(w, 1) < ncomps && @info "$pratio prinicpal ratio achieved with only $(size(w, 1)) components"
    return PCADesign(P, reshape(S, nf, ny, nx), cube, normalize_par_angles(angles))
end

# """
#     pca(cube, angles; ncomps=size(cube, 1), pratio=0.99)

# Use principal component analysis (PCA) to reduce data cube.
# """
# function pca(cube::AbstractArray, angles::AbstractVector, ref::AbstractArray; ncomps = size(cube, 1), pratio = 0.99)
#     # transform cubes
#     rf, ry, rx = size(ref)
#     X = reshape(ref, rf, ry * rx)
#     nf, ny, nx = size(cube)
#     Y = reshape(cube, nf, ny * nx)
    
#     # Fit to our reference, X, and reconstruct
#     P = mvs.fit(mvs.PCA, X; maxoutdim = ncomps, pratio = pratio, mean = 0)
#     S = mvs.reconstruct(P, mvs.transform(P, Y))
    
#     res = cube .- reshape(S, nf, ny, nx)

#     return combine(derotate!(res, angles), method = median)
# end
