struct Framewise{AT} <: ADIAlgorithm
    kernel::AT
    limit
    delta_rot
end

Framewise(alg; limit=Inf, delta_rot=1) = Framewise(alg, limit, delta_rot)

function reconstruct(alg::Framewise, cube::AbstractArray{T,3}; angles, fwhm, r, kwargs...) where T
    pa_threshold = compute_pa_thresh(angles, r, fwhm, alg.delta_rot)
    data = flatten(cube)
    S = similar(data)
    @withprogress name="framewise" begin
        n = 0
        M = size(data, 1)
        @views Threads.@threads for i in axes(data, 1)
            inds = find_angles(angles, i, pa_threshold; limit=alg.limit)
            target = data[i:i, :]
            ref = data[inds, :]
            angs = angles[inds]
            des = fit(alg.kernel, target; kwargs..., ref=ref, angles=angs)
            S[i:i, :] .= reconstruct(des)
            n += 1
            @logprogress n/M
        end
    end
    return expand(S)
end

function reconstruct(alg::Framewise, cube::AnnulusView; ref=cube, angles, fwhm, r=_radius(cube), kwargs...)
    pa_threshold = compute_pa_thresh(angles, r, fwhm, alg.delta_rot)
    data = cube(true) # as view
    S = similar(data)
    @withprogress name="framewise" begin
        n = 0
        M = size(data, 1)
        @views Threads.@threads for i in axes(data, 1)
            inds = find_angles(angles, i, pa_threshold; limit=alg.limit)
            target = data[i:i, :]
            ref = data[inds, :]
            angs = angles[inds]
            des = fit(alg.kernel, target; kwargs..., ref=ref, angles=angs)
            S[i:i, :] .= reconstruct(des)
            n += 1
            @logprogress n/M
        end
    end
    return inverse(cube, S)
end

_radius(cube::AnnulusView) = cube.rmin + (cube.rmax - cube.rmin)/2

function reconstruct(alg::Framewise, cube::MultiAnnulusView; angles, fwhm=cube.width, kwargs...)
    local recons
    anns = eachannulus(cube, true) # as views
    @withprogress name="annulus" begin
        i_ann = 0
        N_ann = length(cube.indices)
        recons = map(anns, cube.radii) do ann, r
            pa_threshold = compute_pa_thresh(angles, r, fwhm, alg.delta_rot)
            S = similar(ann)
            i_angs = 0
            N_angs = length(angles)
            @withprogress name="framewise" begin
                @views Threads.@threads for j in axes(ann, 1)
                    inds = find_angles(angles, j, pa_threshold; limit=alg.limit)
                    target = ann[j:j, :]
                    ref = ann[inds, :]
                    angs = angles[inds]
                    des = fit(alg.kernel, target; kwargs..., ref=ref, angles=angs)
                    S[j:j, :] .= reconstruct(des)
                    i_angs += 1
                    @logprogress i_angs/N_angs
                end
            end
            i_ann += 1
            @logprogress i_ann/N_ann
            return S
        end
    end
    
    return inverse(cube, recons)
end

function reconstruct(alg::Framewise{<:AbstractVector}, cube::MultiAnnulusView; angles, fwhm=cube.width, kwargs...)
    anns = eachannulus(cube, true) # as views
    local recons
    @withprogress name="annulus" begin
        i_ann = 0
        N_ann = length(cube.indices)
        recons = map(anns, cube.radii, alg.kernel) do ann, r, _alg
            pa_threshold = compute_pa_thresh(angles, r, fwhm, alg.delta_rot)
            S = similar(ann)
            @withprogress name="framewise" begin
                i_angs = 0
                N_angs = length(angles)
                @views Threads.@threads for j in axes(ann, 1)
                    inds = find_angles(angles, j, pa_threshold; limit=alg.limit)
                    target = ann[j:j, :]
                    ref = ann[inds, :]
                    angs = angles[inds]
                    des = fit(_alg, target; kwargs..., ref=ref, angles=angs)
                    S[j:j, :] .= reconstruct(des)
                    i_angs += 1
                    @logprogress i_angs/N_angs
                end
            end
            i_ann += 1
            @logprogress i_ann/N_ann
            return S
        end
    end
    return inverse(cube, recons)
end

#######################################

function compute_pa_thresh(angles, r, fwhm, delta_rot)
    pa_threshold = rad2deg(2 * atan(delta_rot * fwhm / 2r))
    mid_range = abs(maximum(angles) - minimum(angles)) / 2
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
