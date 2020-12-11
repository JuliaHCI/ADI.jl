using HCIToolbox: get_annulus_segments

"""
    Median()

Classic PSF subtraction using the median of entire data cube.

# References
1. [Marois et al. 2006](http://adsabs.harvard.edu/abs/2006ApJ...641..556M) Angular Differential Imaging: A Powerful High-Contrast Imaging Technique
"""
struct Median <: ADIAlgorithm end

reconstruct(::Median, cube, angles, cube_ref=cube) = median(cube_ref, dims=1)


################################################################################

# # Optimized median subtraction for a given annulus
# function medsub_annular(cube::AbstractArray{T,3}, angles::AbstractVector; fwhm, ann_size=4, radius_int=0, nframes=4, delta_sep=(0.1, 1), delta_rot=1, kwargs...) where T
#     angles_ = normalize_par_angles(angles)
#     med = median(cube, dims=1)
#     S = cube .- med
#     n, y, x = size(cube)

#     n_annuli = round(Int, (y / 2 - radius_int) / ann_size)

#     # _median_subt_ann_adi
#     ann = 0:n_annuli - 1
#     pa_thresh, inner_r, _ = define_annuli(angles, ann, n_annuli, fwhm, radius_int, ann_size, delta_rot)
#     mask = get_annulus_segments(@view(S[1, :, :]), first(inner_r), ann_size)
#     indices = findall(mask)
#     matrix = S[:, indices]
#     out = similar(matrix)

#     for frame in 1:n
#         if !iszero(pa_thresh)
#             inds_left = find_indices_adi.((angles,), frame, pa_thresh, nframes)
#             @show size(pa_thresh)
#             @show size(inds_left)
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

# function find_indices_adi(angles, frame, pa_thresh, nframes)
#     n = length(angles)

#     index_prev = findfirst(1:frame - 1) do i         
#         abs(angles[frame] - angles[i]) < pa_thresh        
#     end
#     isnothing(index_prev) && (index_prev = frame - 1)

#     index_foll = findfirst(frame:n) do k      
#         abs(angles[k] - angles[frame]) > pa_thresh        
#     end
#     isnothing(index_foll) && (index_foll = n)

#     # nframesi s not None (median annular)
#     window = nframes ÷ 2

#     ind1 = max(index_prev - window, 1)
#     ind2 = index_prev - 1
#     ind3 = index_foll
#     ind4 = min(index_foll + window, n) - 1
#     indices = vcat(ind1:ind2, ind3:ind4)
# end

# """ 
# Function that defines the annuli geometry using the input parameters.
# Returns the parallactic angle threshold, the inner radius and the annulus
# center for each annulus.
# """
# function define_annuli(angles, ann, N, fwhm, radius_int, width, delta_rot)
#     inner_rad = @. ifelse(ann == N - 1, radius_int + (ann * width - 1), radius_int + ann * width)
#     ann_center = @. inner_rad + (width / 2)
#     thresh = @. rad2deg(2 * atan(delta_rot * fwhm / (2 * ann_center)))
#     mid_range = abs(maximum(angles) - minimum(angles)) / 2
#     @. thresh[thresh ≥ 0.9 * mid_range] = 0.9 * mid_range
#     return thresh, inner_rad, ann_center
# end
