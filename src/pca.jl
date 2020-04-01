"""
    pca(cube, angles; ncomps=size(cube, 1), pratio=0.99)

Use principal component analysis (PCA) to reduce data cube.
"""
function pca(cube::AbstractArray, angles::AbstractVector; ncomps = size(cube, 1), pratio = 0.99)
    # transform cube
    nf, ny, nx = size(cube)
    X = reshape(cube, nf, ny * nx)
    
    P = mvs.fit(mvs.PCA, X; maxoutdim = ncomps, pratio = pratio, mean = 0)
    S = mvs.reconstruct(P, mvs.transform(P, X))
    
    res = cube .- reshape(S, nf, ny, nx)

    return combine(derotate!(res, angles), method = median)
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
