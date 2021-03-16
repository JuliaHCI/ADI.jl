
using Distances
using StatsBase: percentile

"""
    LOCI(; dist_threshold=nothing, metric=Cityblock())

Local optimal combination of images (LOCI).

If provided, the frames used for the reference are filtered to only include the frames whose pairwise distances are within the `dist_threshold` percentile. The distances are measured using [Distances.jl](https://github.com/JuliaStats/Distances.jl) and the metric used for measuring the distacnes can be specified with `metric` (by default uses the Manhattan/cityblock metric).

# References

"""
struct LOCI <: ADIAlgorithm
    dist_threshold
    metric::Metric
end

LOCI(; dist_threshold=nothing, metric=Cityblock()) = LOCI(dist_threshold, metric)

function fit(alg::LOCI, data::AbstractMatrix; ref=data, kwargs...)
    coeffs = ref' \ data'
    return LinearDesign(ref, coeffs')
end

function loci_distances_mask(ref::AbstractMatrix, dist_threshold=90, metric=Cityblock())
    distances = pairwise(metric, ref; dims=1)
    thresh = percentile(vec(distances), dist_threshold)
    return @. 0 < distances ≤ thresh
end
loci_distances_mask(ref::AbstractMatrix, ::Nothing, metric) = trues(size(ref)...)

function reconstruct(alg::Framewise{<:LOCI}, cube::AbstractArray{T,3}; angles, fwhm, r, kwargs...) where T
    pa_threshold = compute_pa_thresh(angles, r, fwhm, alg.delta_rot)
    data = flatten(cube)
    S = similar(data)
    dist_mask = loci_distances_mask(data, alg.kernel.dist_threshold, alg.kernel.metric)
    @views Threads.@threads for i in axes(data, 1)
        ang_inds = find_angles(angles, i, pa_threshold; limit=alg.limit)
        inds = ang_inds[dist_mask[i, ang_inds]]
        target = data[i:i, :]
        ref = data[inds, :]
        angs = angles[inds]
        des = fit(alg.kernel, target; kwargs..., ref=ref, angles=angs)
        S[i, :] = reconstruct(des)
    end
    return expand(S)
end

function reconstruct(alg::Framewise{<:LOCI}, cube::AnnulusView; angles, fwhm, r=_radius(cube), kwargs...)
    pa_threshold = compute_pa_thresh(angles, r, fwhm, alg.delta_rot)
    data = cube(true) # as view
    S = similar(data)
    dist_mask = loci_distances_mask(data, alg.kernel.dist_threshold, alg.kernel.metric)
    Threads.@threads for i in axes(S, 1)
        ang_inds = find_angles(angles, i, pa_threshold; limit=alg.limit)
        inds = ang_inds[dist_mask[i, ang_inds]]
        target = data[i:i, :]
        ref = data[inds, :]
        angs = angles[inds]
        des = fit(alg.kernel, target; kwargs..., ref=ref, angles=angs)
        S[i, :] = reconstruct(des)
    end
    return inverse(cube, S)
end

function reconstruct(alg::Framewise{<:LOCI}, cube::MultiAnnulusView; angles, fwhm=cube.width, kwargs...)
    local recons
    anns = eachannulus(cube, true) # as views
    N_ann = length(cube.indices)
    delta_rots = _normalize_deltarot(alg.delta_rot, N_ann)
    @withprogress name="annulus" begin
        i_ann = 0
        recons = map(anns, cube.radii, delta_rots) do ann, r, delta_rot
            pa_threshold = compute_pa_thresh(angles, r, fwhm, delta_rot)
            @debug "PA thresh: $pa_threshold Ann center: $r"
            S = similar(ann)
            dist_mask = loci_distances_mask(ann, alg.kernel.dist_threshold, alg.kernel.metric)
            @views Threads.@threads for j in axes(S, 1)
                ang_inds = find_angles(angles, j, pa_threshold; limit=alg.limit)
                inds = ang_inds[dist_mask[j, ang_inds]]
                target = ann[j:j, :]
                ref = ann[inds, :]
                angs = angles[inds]
                des = fit(alg.kernel, target; kwargs..., ref=ref, angles=angs)
                S[j, :] = reconstruct(des)
            end
            # update progress
            i_ann += 1
            @logprogress i_ann/N_ann
            return S
        end
    end

    return inverse(cube, recons)
end

function reconstruct(alg::Framewise{<:AbstractVector{<:LOCI}}, cube::MultiAnnulusView; angles, fwhm=cube.width, kwargs...)
    anns = eachannulus(cube, true) # as views
    N_ann = length(cube.indices)
    delta_rots = _normalize_deltarot(alg.delta_rot, N_ann)
    local recons
    @withprogress name="annulus" begin
        i_ann = 0
        recons = map(anns, cube.radii, alg.kernel, delta_rots) do ann, r, _alg, delta_rot
            pa_threshold = compute_pa_thresh(angles, r, fwhm, delta_rot)
            S = similar(ann)
            dist_mask = loci_distances_mask(ann, alg.kernel.dist_threshold, alg.kernel.metric)
            @views Threads.@threads for j in axes(S, 1)
                @debug "PA thresh: $pa_threshold Ann center: $r"
                ang_inds = find_angles(angles, j, pa_threshold; limit=alg.limit)
                inds = ang_inds[dist_mask[j, ang_inds]]
                target = ann[j:j, :]
                ref = ann[inds, :]
                angs = angles[inds]
                des = fit(_alg, target; kwargs..., ref=ref, angles=angs)
                S[j, :] = reconstruct(des)
            end
            # update progress
            i_ann += 1
            @logprogress i_ann/N_ann
            return S
        end
    end
    return inverse(cube, recons)
end

