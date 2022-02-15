
using Distances
using StatsBase: quantile

"""
    LOCI(; dist_threshold=nothing, metric=Cityblock())

Local optimal combination of images (LOCI).

If provided, the frames used for the reference are filtered to only include the frames whose pairwise distances are within the `dist_threshold` quantile. The distances are measured using [Distances.jl](https://github.com/JuliaStats/Distances.jl) and the metric used for measuring the distacnes can be specified with `metric` (by default uses the Manhattan/cityblock metric).

!!! note "Traditional LOCI"
    Unless [`LOCI`](@ref) is applied [`Framewise`](@ref), no distance thresholding will occur. In order to recreate the traditional LOCI algorithm, consider constructing an algorithm like
    ```julia
    # only use the most-similar 90th-percentile frames
    alg = Framewise(LOCI(dist_threshold=0.90))
    ```

# References
1. [Lafreniere et al. (2007)](http://adsabs.harvard.edu/abs/2007ApJ...660..770L) "A New Algorithm for Point-Spread Function Subtraction in High-Contrast Imaging: A Demonstration with Angular Differential Imaging"
"""
@concrete struct LOCI{M<:Metric} <: ADIAlgorithm
    dist_threshold
    metric::M
end

LOCI(; dist_threshold=nothing, metric::Metric=Cityblock()) = LOCI(dist_threshold, metric)

function fit(alg::LOCI, data::AbstractMatrix; ref=data, kwargs...)
    coeffs = ref \ data
    return LinearDesign(ref, coeffs)
end

function loci_distances_mask(ref::AbstractMatrix, dist_threshold=0.90, metric=Cityblock())
    distances = pairwise(metric, ref; dims=2)
    thresh = quantile(vec(distances), dist_threshold)
    return @. 0 < distances ≤ thresh
end
loci_distances_mask(ref::AbstractMatrix, ::Nothing, metric) = Fill(true, size(ref))

function reconstruct(alg::Framewise{<:LOCI}, cube::AbstractArray{T,3}; angles, kwargs...) where T
    pa_threshold = compute_pa_thresh(angles, alg.delta_rot; kwargs...)
    data = flatten(cube)
    S = similar(data)
    dist_mask = loci_distances_mask(data, alg.kernel.dist_threshold, alg.kernel.metric)
    @views Threads.@threads for i in axes(S, 2)
        ang_inds = find_angles(angles, i, pa_threshold; limit=alg.limit)
        inds = ang_inds[dist_mask[i, ang_inds]]
        target = data[:, i:i]
        ref = data[:, inds]
        angs = angles[inds]
        des = fit(alg.kernel, target; kwargs..., ref=ref, angles=angs)
        S[:, i] = reconstruct(des)
    end
    return expand(S)
end

function reconstruct(alg::Framewise{<:LOCI}, cube::AnnulusView; angles, r=_radius(cube), kwargs...)
    pa_threshold = compute_pa_thresh(angles, alg.delta_rot; r=r, kwargs...)
    data = cube(true) # as view
    S = similar(data)
    dist_mask = loci_distances_mask(data, alg.kernel.dist_threshold, alg.kernel.metric)
    Threads.@threads for i in axes(S, 2)
        ang_inds = find_angles(angles, i, pa_threshold; limit=alg.limit)
        inds = ang_inds[dist_mask[ang_inds, i]]
        target = data[:, i:i]
        ref = data[:, inds]
        angs = angles[inds]
        des = fit(alg.kernel, target; kwargs..., ref=ref, angles=angs)
        S[:, i] = reconstruct(des)
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
            pa_threshold = compute_pa_thresh(angles, delta_rot; fwhm=fwhm, r=r)
            @debug "PA thresh: $pa_threshold Ann center: $r"
            S = similar(ann)
            dist_mask = loci_distances_mask(ann, alg.kernel.dist_threshold, alg.kernel.metric)
            @views Threads.@threads for j in axes(S, 2)
                ang_inds = find_angles(angles, j, pa_threshold; limit=alg.limit)
                inds = ang_inds[dist_mask[ang_inds, j]]
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

function reconstruct(alg::Framewise{<:AbstractVector{<:LOCI}}, cube::MultiAnnulusView; angles, fwhm=cube.width, kwargs...)
    anns = eachannulus(cube, true) # as views
    N_ann = length(cube.indices)
    delta_rots = _normalize_deltarot(alg.delta_rot, N_ann)
    local recons
    @withprogress name="annulus" begin
        i_ann = 0
        recons = map(anns, cube.radii, alg.kernel, delta_rots) do ann, r, _alg, delta_rot
            pa_threshold = compute_pa_thresh(angles, delta_rot; r=r, fwhm=fwhm)
            S = similar(ann)
            dist_mask = loci_distances_mask(ann, _alg.dist_threshold, _alg.metric)
            @views Threads.@threads for j in axes(S, 2)
                @debug "PA thresh: $pa_threshold Ann center: $r"
                ang_inds = find_angles(angles, j, pa_threshold; limit=alg.limit)
                inds = ang_inds[dist_mask[ang_inds, j]]
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

