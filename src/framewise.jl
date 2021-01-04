struct Framewise{AT,OT} <: ADIAlgorithm
    kernel::AT
    opts::OT
end

Framewise(alg::ADIAlgorithm; options...) = Framewise(alg, options)

function reconstruct(alg::Framewise, cube::AbstractArray{T,3}; angles, fwhm, r, kwargs...) where T
    pa_threshold = compute_pa_thresh(angles, r, fwhm, delta_rot)
    data = flatten(cube)
    S = similar(data)
    N = length(angles)
    p = Progress(N; desc="framewise ")
    Threads.@threads for i in axes(data, 1)
        inds = find_angles(angles, i, pa_threshold; alg.limit)
        target = @view data[i, :]
        ref = @view data[inds, :]
        angs = @view angles[inds]
        S[i, :] .= reconstruct(alg.kernel, target; kwargs..., ref=ref, angles=angs)
        next!(p)
    end
    return expand(S)
end

function reconstruct(alg::Framewise, cube::AnnulusView; ref=cube, angles, fwhm, r=(cube.rmax - cube.rmin)/2, kwargs...)
    pa_threshold = compute_pa_thresh(angles, r, fwhm, delta_rot)
    data = cube(true)
    S = similar(data)
    p = Progress(N; desc="framewise ")
    Threads.@threads for i in axes(data, 1)
        inds = find_angles(angles, i, pa_threshold; alg.limit)
        target = @view data[i, :]
        ref = @view data[inds, :]
        angs = @view angles[inds]
        S[i, :] .= reconstruct(alg.kernel, target; kwargs..., ref=ref, angles=angs)
        next!(p)
    end
    return inverse(cube, S)
end

# function reconstruct(alg::Framewise, cube::MultiAnnulusView; ref=cube, angles, fwhm=cube.width, kwargs...)
#     pa_threshold = compute_pa_thresh(angles, r, fwhm, delta_rot)
#     Threads.@threads for i in axes(data, 1)
#         inds = find_angles(angles, i, pa_threshold; alg.limit)
#     end
# end

#######################################

function compute_pa_thresh(angles, r, fwhm, delta_rot)
    pa_threshold = rad2deg(2atan(delta_rot * fwhm / (2r)))
    mid_range = abs(-(extrema(angles)...)) / 2
    k = mid_range * 0.9
    if pa_threshold ≥ k
        @info "pa threshold $pa_threshold too large, will be set to $k"
        pa_threshold = k
    end
    return pa_threshold
end

function find_angles(angles, idx, thresh; limit=Inf)
    fidx, lidx = firstindex(angles), lastindex(angles)
    p, n = fidx, idx
    for i in fidx:idx-1
        if abs(angles[idx] - angles[i]) < thresh
            p = i
            break
        else
            p += 1
        end
    end
    for k in idx:lidx
        if abs(angles[k] - angles[idx]) > thresh
            n = k
            break
        else
            n += 1
        end
    end
    if isfinite(limit)
        window = limit ÷ 2
        first_half = max(p - window, fidx):p - 1
        last_half = n:min(n + window, lidx) - 1
    else
        first_half = fidx:p - 1
        last_half = n:lidx
    end
    return vcat(first_half, last_half)
end