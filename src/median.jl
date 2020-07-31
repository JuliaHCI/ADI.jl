"""
    Median()

Classic PSF subtraction using the median of entire data cube.

# References
* [Marois et al. 2006](https://arxiv.org/abs/astro-ph/0512335) Angular Differential Imaging: A Powerful High-Contrast Imaging Technique
"""
struct Median <: ADIAlgorithm end

function (::Median)(cube, angles, ref=cube, args...; kwargs...)
    angles_ = normalize_par_angles(angles)
    med = median(ref, dims=1)
    S = cube .- med
    return collapse!(S, angles_; kwargs...)
end

# # Optimized median subtraction for a given annulus
# function medsub_annular(cube::AbstractArray{T,3}, angles::AbstractVector; fwhm, ann_size=4, radius_int=0, nframes=4, delta_sep=(0.1, 1), delta_rot=1, kwargs...)
#     angles_ = normalize_par_angles(angles)
#     med = median(cube, dims=1)
#     S = cube .- med
#     n, y, x = size(cube)

#     n_annuli = round(Int, (y / 2 - radius_int) / ann_size)

#     # _median_subt_ann_adi
#     pa_thresh, inner_r, _ = define_annuli(angles, ann, N, fwhm, radius_int, width, delta_rot, 1, false)

#     mask = get_annulus_segments(view(S[1, :, :]), inner_r, width) |> first
#     indices = findall(mask)
#     matrix = S[:, indices]
#     out = similar(matrix)

#     for frame in axes(matrix, 1)
#         if !iszero(pa_thresh)
#             inds_left = find_indices_adi(angles, frame, pa_thresh, nframes)
#             matrix_disc = @view matrix[inds_left, :, :]
#         else
#             matrix_disc = matrix
#         end
#         ref_psf_opt = median(matrix_disc, dims=1)
#         curr_frame = @view matrix[frame, :, :]
#         @. out[frame, :, :] = curr_frame - ref_psf_opt
#     end
#     # return out, yy, xx, pa_thresh
#     # end _median_subt_ann_adi

#     for  ann in 1:n_annuli
#         out[:, yy[ann], xx[ann]] .= out[ann, :, :]
#     end

# end

# function find_indices_adi(angles, frame, pa_thresh, nframe)
#     n = length(angles)
#     index_prev = firstindex(angles)
#     index_foll = frame
#     for i in 1:frame
#         if abs(angles[frame] - angles[i]) < pa_thresh
#             index_prev = i
#             break
#         else
#             index_prev += 1
#         end
#     end

#     for k in frame:n
#         if  abs(angles[k]  - angles[frame]) > pa_thresh
#             index_foll = k
#             break
#         else
#             index_foll += 1
#         end
#     end

#     # nframesi s not None (median annular)
#     window = nframes ÷ 2

#     ind1 = max(index_prev - window, 1)
#     ind2 = index_prev
#     ind3 = index_foll
#     ind4 = min(index_foll + window, n)
#     indices = vcat(ind1:ind2, ind3:ind4)
# end
