using Statistics

"""
    stim(residuals, angles)

Calculate the standardized trajectory intensity mean (STIM) map. The inputs are a cube of residuals and the corresponding parallactic angles.

This method seeks to improve upon the typical student-t S/N tests ([`snr`](@ref), [`significance`](@ref)) by calculating statistics in the temporal domain instead of the spatial domain. This is why the full residual cube is required rather than a reduced frame. 

In particular, the STIM map is robust to detections with multiple objects or extended sources within the same annuli, which results in very high noise estimates using spatial methods. The STIM map also performs better at small angular separations, since the temporal domain has no limitations from limited resolution elements.

*Pairet et al. 2019* derives a detection threshold of `τ ≈ 0.5` for the STIM map. The detection threshold can be calculated in a similar manner using [`stim_threshold`](@ref).

# Examples

```julia
julia> cube, angles = # load data

julia> L = reconstruct(PCA(10), cube, angles);

julia> S = cube .- L;

julia> stimmap = stim(S, angles);
```

# References
* [Pairet et al. 2019](http://adsabs.harvard.edu/abs/2019MNRAS.487.2262P) "STIM map: detection map for exoplanets imaging beyond asymptotic Gaussian residual speckle noise"

# See Also
[`stim_threshold`](@ref)
"""
function stim(residuals::AbstractArray{T,3}, angles) where T
    return stimmap(derotate(residuals, angles))
end

"""
    stim_threshold([stimmap, ] residuals, angles)

Calculate the detection threshold for the standardized trajectory intensity mean (STIM) map. This method uses the same residual cube as [`stim`](@ref) but adds an additional step of estimating the residual noise by derotating the residuals with the *opposite* parallactic angles.

If the STIM map has already been calculated, it can be passed in, otherwise it will be calculated in addition to the noise map. Note this will not return the STIM map, only the threshold.

The threshold is derived in section 5.1 of *Pairet et al. 2019* as the ratio of the number of pixels above the approximated noise map.

# References
* [Pairet et al. 2019](http://adsabs.harvard.edu/abs/2019MNRAS.487.2262P) "STIM map: detection map for exoplanets imaging beyond asymptotic Gaussian residual speckle noise"

# See Also
[`stim`](@ref)
"""
function stim_threshold(stimmap, residuals, angles)
    # estimate noise map by getting STIM map of cube
    # rotated with opposite angles
    d_opp = Metrics.stimmap(derotate(residuals, -angles))
    # return ratio of values above the noise
    n_ϵ = count(stimmap .> d_opp)
    n = length(stimmap)
    return n_ϵ / n
end

function stim_threshold(residuals, angles)
    # calculate stimmap first
    d = stimmap(derotate(residuals, angles))
    return stim_threshold(d, residuals, angles)
end

"""
    Metrics.stimmap(cube)

Calculates the STIM map of a derotated cube. This is roughly equivalent to the temporal mean divided by the temporal standard deviation.
"""
function stimmap(cube)
    μ = mean(cube, dims=1)
    σ = std(cube, dims=1, mean=μ)
    d = zero(μ)
    @. @views d[!iszero(σ)] = μ[!iszero(σ)] / σ[!iszero(σ)]
    return dropdims(d, dims=1)
end
