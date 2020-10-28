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
       throughput

include("contrast.jl")
include("snr.jl")

end
