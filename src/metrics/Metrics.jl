"""
    ADI.Metrics

This module provides code for analyzing the results from ADI in a way that is interpretable statistically. Some of the key functionalities are signal-to-noise, significance, the receiver operating characteristic, and the contrast curve.
"""
module Metrics

export detectionmap,
       snr,
       significance,
       noise,
       contrast_curve,
       throughput,
       stim,
       stim_threshold

include("contrast.jl")
include("snr.jl")
include("stim.jl")

end
