using ProgressLogging

# struct PairetDesign{D<:ADIDesign,T<:AbstractArray,V<:AbstractVector} <: ADIDesign{T, V}
#     des::D{T, V}
# end

# function Base.show(io::IO, d::PairetDesign{D,T}) where {D,T}
#     n = mvs.indim(d.des.pca)
#     p = mvs.outdim(d.des.pca)
#     pr = mvs.principalratio(d.des.pca)
#     print("PairetDesign{$(D.name){$T}}(ncomps=$p, D=$n, pratio=$pr)")
# end


pairet(cube::AbstractArray, angles::AbstractVector; ncomps, threshold=0) = 
    pairet(:pca, cube, angles; ncomps=ncomps, threshold=threshold)

pairet(s::Symbol, cube, angles; ncomps, threshold=0) = 
    pairet(Val(s), cube, angles; ncomps=ncomps, threshold=threshold)


function pairet(::Union{Val{:pca}, Val{:PCA}}, cube::AbstractArray, angles::AbstractVector; ncomps, threshold=0)
    des = pca(cube, angles, ncomps=1, mean=nothing)
    @progress for n in 1:ncomps
        red = reduce(des, cube, des._angs)
        R = cube .- _pairet_theta(red, des._angs, threshold)
        des = pca(R, des._angs, ncomps=n, mean=nothing)
        des.S .= _pairet_cube_est(des, cube)
    end
    des._cube .= cube
    return des
end

# takes a frame, expands it into a cube, rotates it clockwise by angles, 
# and min-clips at threshold
function _pairet_theta(frame, angles, threshold)
    N = length(angles)
    cube = similar(frame, N, size(frame)...)
    for idx in axes(cube, 1)
        @. cube[idx, :, :] = frame
    end
    @. cube[cube ≤ threshold] = 0
    return derotate!(cube, -angles)
end

function _pairet_cube_est(des, cube)
    m = mean(des._cube, dims=1)
    target = flatten(cube .- m)
    zvectors = mvs.transform(des.pca, flatten(cube))

    best_cube = similar(target)
    for i in axes(best_cube, 1)
        v11 = sum(zvectors .* target[i:i, :], dims=2)
        best_cube[i, :] .= vec(v11' * zvectors)
    end
    return m .+ expand(best_cube)
end

# pairet(cube::AbstractArray, angles::AbstractVector; ncomps) = pairet(pca, cube, angles; ncomps=ncomps)
