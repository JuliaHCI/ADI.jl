"""
    pca(cube, angles; ncomps=size(cube, 1), pratio=0.99)

Use principal component analysis (PCA) to reduce data cube.

# Examples

```jldoctest
julia> cube = ones(30, 100, 100);

julia> d = design(ADI.PCA, cube);

julia> d.info
Dict{Symbol,Any} with 2 entries:
  :ncomps => 30
  :pratio => 1.0

julia> d = design(ADI.PCA, cube, ncomps=5);

julia> d.info
Dict{Symbol,Any} with 2 entries:
  :ncomps => 5
  :pratio => 1.0
```
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
