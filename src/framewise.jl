@concrete struct Framewise <: ADIAlgorithm
    kernel
    limit
    delta_rot
end

"""
    Framewise(alg; limit=Inf, delta_rot=nothing)
    Framewise(algs::AbstractVector; limit=Inf, delta_rot=nothing)

Wrap an algorithm such that the underlying data will be processed frame by frame. For each frame a reference library is created from the data. This reference can be filtered by rejecting frames which have not rotated a sufficient parallactic angle. `delta_rot` sets the required arc length for rotation in units of the FWHM. If `delta_rot` is `nothing` there will be now temporal filtering. The number of frames retained can be specified with `limit`, e.g. the 4 closest frames in time with the target frame.

The following keyword arguments must be provided to [`reconstruct`](@ref) or [`subtract`](@ref)
* `angles` - the measured parallactic angles for each frame

If `delta_rot` is provided, the following additional keyword arguments must be provided to [`reconstruct`](@ref) or [`subtract`](@ref).
* `fwhm` - the FWHM of the instrument in pixels. Will be set to the width of a [`MultiAnnulusView`](@ref)
* `r` - The radius of the arc to calculate the parallactic angle threshold. Will be set automatically if using [`AnnulusView`](@ref) or [`MultiAnnulusView`](@ref).

In addition, `Framewise` versions of algorithms do not implement [`ADI.fit`](@ref) and do not currently support RDI.

## Annular reduction

In the case of reducing a [`MultiAnnulusView`](@ref), a vector of algorithms can be used, each one corresponding to each annulus of the view. In this case, too, `delta_rot` can be given as a vector or as a tuple. If it is given as a tuple, `delta_rot` will increase linearly from the first value to the last value across each annulus.

# Examples

```julia
julia> cube, angles = # load data

julia> alg = PCA(10) # the algorithm to use on each reference

julia> res = Framewise(alg)(cube, angles);

julia> mav = MultiAnnulusView(cube, 5; inner=5);

julia> res_ann = Framewise(alg, delta_rot=(0.1, 1))(mav, angles);
```
"""
Framewise(alg; limit=Inf, delta_rot=nothing) = Framewise(alg, limit, delta_rot)

function reconstruct(alg::Framewise, cube::AbstractArray{T,3}; angles, kwargs...) where T
    pa_threshold = compute_pa_thresh(angles, alg.delta_rot; kwargs...)
    data = flatten(cube)
    S = similar(data)
    @views Threads.@threads for i in axes(data, 2)
        inds = find_angles(angles, i, pa_threshold; limit=alg.limit)
        target = data[:, i:i]
        ref = data[:, inds]
        angs = angles[inds]
        des = fit(alg.kernel, target; kwargs..., ref=ref, angles=angs)
        S[:, i] = reconstruct(des)
    end
    return expand(S)
end

function reconstruct(alg::Framewise, cube::AnnulusView; angles, r=_radius(cube), kwargs...)
    pa_threshold = compute_pa_thresh(angles, alg.delta_rot; r=r, kwargs...)
    data = cube(true) # as view
    S = similar(data)
    Threads.@threads for i in axes(S, 2)
        inds = find_angles(angles, i, pa_threshold; limit=alg.limit)
        target = data[:, i:i]
        ref = data[:, inds]
        angs = angles[inds]
        des = fit(alg.kernel, target; kwargs..., ref=ref, angles=angs)
        S[:, i] = reconstruct(des)
    end
    return inverse(cube, S)
end

_radius(cube::AnnulusView) = cube.rmin + (cube.rmax - cube.rmin)/2

function reconstruct(alg::Framewise, cube::MultiAnnulusView; angles, fwhm=cube.width, kwargs...)
    local recons
    anns = eachannulus(cube, true) # as views
    N_ann = length(cube.indices)
    delta_rots = _normalize_deltarot(alg.delta_rot, N_ann)
    @withprogress name="annulus" begin
        i_ann = 0
        recons = map(anns, cube.radii, delta_rots) do ann, r, delta_rot
            pa_threshold = compute_pa_thresh(angles, delta_rot; r=r, fwhm=fwhm)
            @debug "PA thresh: $pa_threshold Ann center: $r"
            S = similar(ann)
            @views Threads.@threads for j in axes(S, 2)
                inds = find_angles(angles, j, pa_threshold; limit=alg.limit)
                target = ann[:, j:j]
                ref = ann[:, inds]
                angs = angles[inds]
                des = fit(alg.kernel, target; kwargs..., ref=ref, angles=angs)
                S[:, j] = reconstruct(des)
            end
            # update progress
            i_ann += 1
            @logprogress i_ann/N_ann
            return S
        end
    end

    return inverse(cube, recons)
end

_normalize_deltarot(delta_rot, N) = Fill(delta_rot, N)
_normalize_deltarot(delta_rot::Tuple, N) = range(delta_rot..., length=N)
_normalize_deltarot(delta_rot::AbstractVector, N) = delta_rot

function reconstruct(alg::Framewise{<:AbstractVector}, cube::MultiAnnulusView; angles, fwhm=cube.width, kwargs...)
    anns = eachannulus(cube, true) # as views
    N_ann = length(cube.indices)
    delta_rots = _normalize_deltarot(alg.delta_rot, N_ann)
    local recons
    @withprogress name="annulus" begin
        i_ann = 0
        recons = map(anns, cube.radii, alg.kernel, delta_rots) do ann, r, _alg, delta_rot
            pa_threshold = compute_pa_thresh(angles, delta_rot; fwhm=fwhm, r=r)
            S = similar(ann)
            @views Threads.@threads for j in axes(S, 2)
                @debug "PA thresh: $pa_threshold Ann center: $r"
                inds = find_angles(angles, j, pa_threshold; limit=alg.limit)
                target = ann[:, j:j]
                ref = ann[:, inds]
                angs = angles[inds]
                des = fit(_alg, target; kwargs..., ref=ref, angles=angs)
                S[:, j] = reconstruct(des)
            end
            # update progress
            i_ann += 1
            @logprogress i_ann/N_ann
            return S
        end
    end
    return inverse(cube, recons)
end

#######################################

compute_pa_thresh(angles, delta_rot::Nothing; kwargs...) = nothing

function compute_pa_thresh(angles, delta_rot; r, fwhm, kwargs...)
    pa_threshold = 2 * atand(delta_rot * fwhm,  2r)
    mid_range = abs(maximum(angles) - minimum(angles)) / 2
    k = mid_range * 0.9
    if pa_threshold ≥ k
        @info "pa threshold $pa_threshold too large, will be set to $k"
        pa_threshold = k
    end
    return pa_threshold
end

function find_angles(angles, idx, thresh::Nothing; limit=Inf)
    isfinite(limit) || return vcat(1:idx - 1, idx + 1:length(angles))
    window = limit ÷ 2
    first_half = max(idx - 1 - window, 1):idx - 1
    last_half = idx + 1:min(idx + 1 + window, length(angles))
    return vcat(first_half, last_half)
end
function find_angles(angles, idx, thresh; limit=Inf)
    iszero(thresh) && return vcat(1:idx - 1, idx + 1:length(angles))
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
