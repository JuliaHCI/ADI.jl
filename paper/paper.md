---
title: 'ADI.jl: A Julia Package for High-Contrast Imaging'
tags:
  - Julia
  - astronomy
  - high-contrast imaging
  - direct imaging
  - image processing
authors:
  - name: Miles Lucas
    orcid: 0000-0001-6341-310X
    affiliation: 1
  - name: Michael Bottom
    orcid: 0000-0003-1341-5531
    affiliation: 1
affiliations:
 - name: Institute for Astronomy, University of Hawai'i
   index: 1
date: 11/02/2020
bibliography: paper.bib
---

# Summary

High-contrast imaging (HCI) is a powerful technique for discovering and characterizing exoplanets. Being able to probe the architecture, formation, and atmospheres of planets *directly* is necessary for advancing companion formation and evolution theory [@bowler_imaging_2016]. The process required to image an exoplanet is daunting, however, due to the brightness and proximity of their host stars. One part of making such a difficult detection is the image processing of HCI data to remove systematic signal and attenuate noise. The size of HCI data and complexity of processing algorithms requires an efficient numerical framework while simultaneously offering modularity for rapid experimentation and discovery.

# Statement of need

`ADI.jl` is a Julia framework for the post-processing of HCI data. The high-level composability yet low-level efficiency of Julia [@Julia-2017] is well suited for this task. `ADI.jl` organizes common post-processing patterns into separate functions that are overloaded by new algorithms. In addition, the techniques for applying these algorithms, such as full-frame reduction versus annular reduction, are completely decoupled, functionally, from the algorithms. This creates an API that is highly modular and benefits both observers and algorithm designers. Observers can quickly try different algorithms and application methods for the rapid discovery of exoplanets. Julia's multiple dispatch [@2018arXiv180803370B] and rich language interoperability  means developers have to write less code (in their language of choice) to implement the `ADI.jl` API for a new algorithm.

# A generic framework for high-contrast imaging

The data in high-contrast imaging is complicated and requires significant pre-processing and post-processing to enable the discovery and characterization of exoplanets. Post-processing for HCI comprises estimating the systematic noise and removing it. Each HCI algorithm mainly differs in how this systematic signal is estimated. Fortunately, though, the applications of these algorithms follow a common process that forms the generic framework for `ADI.jl`.

For a brief primer, angular differential imaging (ADI) is an observational technique that utilizes the Earth's rotation throughout a night of observing [@liu:2004; @marois_angular_2006]. Normally telescopes have optics to counter this rotation, but disabling these optics and taking images throughout the night will give a sequence of frames with rotatinng astrophysical signal. The systematic signal, though, is mainly affected by the telescope optics, so the systematics do not appear to rotate like the sky. The sequence of frames are aligned and stacked together into a data cube. ADI algorithms exploit the difference in rotation to approximate and subtract the systematics. The remaining sequence of residual images can then be reoriented and combined to attenuate the noise and reveal potential companions.

In `ADI.jl` this process is implemented in a functional framework. Each algorithm implements a `reconstruct` method, which approximates the systematic signal from a *target* data cube or using a separate *reference* data cube (referred to as reference differential imaging, or RDI). The algorithm is a first-class function, which means the options specific to the algorithm can be separated from the arguments for the usage of the algorithm. For example, the common reduction algorithm of principal component analysis (PCA) [@soummer_detection_2012] requires specifying the number of components to fit to the data. This parameter is set within the `PCA` algorithm which is then used generically with the `ADI.jl` framework. In other words, the specification of the algorithm is completely separated from the application of the algorithm.

The applications of these algorithms can vary from unsupervised end-to-end pipelines to extensively hand-tuned reductions seeking the best signal-to-noise ratio possible. Recognizing the different applications of the algorithms guides our design of the `ADI.jl` API. Algorithms can be used for full reductions in a single line of code, the speckle estimate can be obtained using `reconstruct` and methods like PCA, which have a compact linear form, have a `decompose` method for the decomposition of the speckle estimate into a basis and weights. So far this level of granularity has proven appropriate. For example, the fixed-point greedy disk subtraction algorithm (GreeDS) [@pairet_reference-less_2019; @pairet_mayonnaise_2020] iteratively forms a linear basis and applies a non-linear correction. The details of this algorithm are outlined using PCA as the method for obtaining the basis, but there is no reason it cannot be applied with non-negative matrix factorization (NMF) [@ren_non-negative_2018]. Therefore the `GreeDS` algorithm utilizes the `decompose` method *generically*, giving users flexibility without forcing developers to write out all possible combinations of methods.

# The unreasonable effectiveness of multiple-dispatch

In addition to the various ADI algorithms there are different techniques for how the data is presented. Full-frame reduction uses the entire sequence of images for fitting the PSF. Alternative techniques, such as dividing the data spatially into annuli or sectors, temporal filtering, and computational methods such as out-of-core processing or using GPUs all offer an extra dimension of experimentation for observers and an extra dimension of implementation details for algorithm designers. For every new algorithm or technique the combination of inputs increases exponentially. This exponential growth is detrimental for developers because each new experiment requires copying or re-implementing large quantities of code, or otherwise abandoning use-cases not immediately applicable. `ADI.jl` overcomes this burden uniquely through the use of multiple-dispatch.

Multiple-dispatch allows functions to dispatch on the types of all input arguments. This is in contrast to common single-dispatch languages, like C and Python, which can only organize their algorithms using inheritance or function overloading^[function overloading is distinct from multiple-dispatch because function overloading only dispatches on the overloaded caller, not on any of the other types. This is roughly the same reason why mix-ins or monkey patching are not multiple-dispatch.]. For a given set of input argument types, single-disptach provides a linear number of methods for a single function, while multiple-dispatch provides an exponential number of methods. Applying this to HCI, if the algorithms obey a type structure and the delivery methods obey a type structure, only a linear number of methods have to be written in Julia, while an exponential number of methods would have to be written in a single-dispatch language like Python.

What this means for users is that mixing and matching different pieces of algorithms and delivery methods just works. For example, if `ADI.jl` knows *how* to do annular reduction regardless of the algorithm, then multiple-dispatch enables any algorithm that appropriately implements the API described in the previous section to be applied, to the annuli or sectors. This is very useful for developers because the core essence of the algorithm is all that needs to be implemented, without having to worry about delivery methods or other interactions. In addition, Julia's foreign function interface (FFI) and language interoperability allows using code written in Python or C for example, directly within `ADI.jl`. This culminates in an HCI framework that can be used with existing packages and tools while offering a strong base for developing and extending algorithms.

# Comparisons with existing software

High-contrast imaging, as a field, predominantly utilizes Python for data reduction. We break down some of the necessary computations into *pre-processing*, which includes raw calibration of data, the centering and stacking of the data cube, bad-pixel removal, etc., *post-processing*, which includes the PSF approximation and subtraction, *posterior analysis* which includes analyzing the results of post-processing like detection statistics and throughput calibration, and finally *forward modeling* which includes various statistical models for companions and disks for use with post-processing algorithms. `ADI.jl` primarily focuses on post-processing and posterior analysis.

Some notable libraries for HCI tasks include the Vortex Imaging Pipeline (`VIP`) [@gomez_gonzalez_vip_2016], `pyKLIP` [@2015ascl.soft06001W], and `PynPoint` [@pynpoint:2019]. A table of the feature sets of these packages alongside `ADI.jl` is presented in \autoref{tab:features}. In particular,`VIP`has served as a useful source of information regarding HCI image-processing as well as detailed implementations of common ADI algorithms. This has been indispensable in the development of `ADI.jl`, although this package is not a direct translation.

`VIP` offers the most diversity in algorithms and their applications, but at the cost of flexibility. The implementations of the algorithms suffer from long lists of arguments, large blocks of parsing the input arguments and coercing them into usable types or manually dispatching. Because of this, algorithms with similar techniques, such as annular processing using the median vs. annular processing with PCA, are not guaranteed to be applied the same way and any new algorithms must repeat code and add complexity in order to use the annular processing technique.

Table 1: Comparison of features across different HCI frameworks \label{tab:features}

 Framework | Pre-processing | Algorithms | Techniques | Posterior analysis | Forward modeling
-|-|-|-|-|-
ADI.jl | no | median, PCA, NMF, fixed-point GreeDS | Full-frame ADI/RDI, SDI (experimental) | detection maps, contrast curve | no
VIP | yes | median, LOCI, PCA, NMF, LLSG, ANDROMEDA | Full-frame ADI/RDI, SDI, annular ADI/RDI for some algorithms | detection maps, blob detection, STIM, ROC, contrast curve | NegFC
pyKLIP | no | PCA | Full-frame ADI/RDI, SDI, annular ADI/RDI | detection maps, blob detection, contrast curve | KLIP-FM, Planet Evidence, matched filter (FMMF), spectrum fitting, DiskFM
PynPoint | yes | median, PCA | Full-frame ADI/RDI, SDI | detection maps, contrast curve | no

# Acknowledgments

We acknowledge the `VIP` team for creating a package that has been an indispensable resource in the creation of `ADI.jl`. We also acknowledge the ongoing research using `ADI.jl` for forward modeling of exoplanets with the package `Firefly.jl`.

# References
