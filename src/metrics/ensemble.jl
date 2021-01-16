
"""
    slimmap(residuals, angles; N)
    slimmap(stim_maps; N)

Calculate the STIM largest intensity mask (SLIMask) map. This is an ensemble method which averages the STIM maps for multiple residual cubes. In addition to computing this average, for each STIM map all but the brightest `N` pixels will be masked out and eventaully averaged to create the SLIMask.

`N` should represent a cutoff for the number of expected companions. For example, if the FWHM of the companion signal is 5 pixels, then the area under the fwhm is ~20 pixels. If I want to probe the brightest 4 potential companions, I would mask all but the `N = 20 * 4 = 80` brightest pixels.

Both the average STIM map and the SLIMask will be returned, and the two can be multiplied together to produce the SLIMask-STIM map. This achieves two things: first, the masking makes most of the pixels 0, providing better visual contrast in the map, and second, by averaging the mask, pixels which are not consistently in the brightest `N` for each STIM map will have lower probabilities in the corresponding SLIMask-STIM map.

# Examples

This example recreates the analysis shown in Pairet, B. (2020) where the SLIM map is computed with the ensemble of residual cubes produced by increasing ranks of PCA subtraction.

```julia
julia> cube, angles = # load data

julia> algs = PCA.(5:25);

julia> residual_cubes = subtract(algs, cube);

julia> stim_av, slimmask = slimmap(residual_cubes, angles; N=100);

julia> slim_prob_map = stim_av .* slimmask;
```

# References

1. [Pairet, B. 2020](https://dial.uclouvain.be/pr/boreal/object/boreal:240621) "Signal processing methods for high-contrast observations of planetary systems"

# See also

[`stimmap`](@ref), [`stim`](@ref)
"""
function slimmap(cubes::AbstractVector{CT}, angles; N) where {T,CT<:AbstractArray{T,3}}
    # simulatenous calculate mask and mean map
    stim_av = similar(first(cubes), size(first(cubes))[2:3]...)
    mask_av = similar(stim_av)
    η = 1 / length(cubes)
    Threads.@threads for residual in cubes
        stim_map = stimmap(residual, angles)
        stim_av .+= stim_map
        # get Nth brightest pixel as mask threshold
        thresh = partialsort(vec(stim_map), N; rev=true)
        @. mask_av[stim_map ≥ thresh] += η
    end
    # normalize to get mean
    stim_av .*= η
    return stim_av, mask_av
end

function slimmap(stimmaps::AbstractVector{ST}; N) where {ST<:AbstractMatrix}
    # simulatenous calculate mask and mean map
    stim_av = similar(first(stimmaps))
    mask_av = zero(stim_av)
    η = 1 / length(stimmaps)
    Threads.@threads for stim_map in stimmaps
        stim_av .+= stim_map
        # get Nth brightest pixel as mask threshold
        thresh = partialsort(vec(stim_map), N; rev=true)
        @. mask_av[stim_map ≥ thresh] += η
    end
    # normalize to get mean
    stim_av .*= η
    return stim_av, mask_av
end
