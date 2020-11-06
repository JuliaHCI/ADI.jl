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

High-contrast imaging (HCI) is a powerful technique for discovering and characterizing exoplanets. Being able to probe the architecture, formation, and atmospheres of planets *directly* is necessary for advancing companion formation and evolution theory [@bowler_imaging_2016]. The process required to image an exoplanet is daunting, however, due to the brightness and proximity of their host stars. One part of making such a difficult detection is the image processing of HCI data to remove systematic signal and attenuate noise. The size of HCI data and complexity of processing algorithms requires an efficient numerical framework that simultaneously offers modularity for rapid experimentation and discovery.

The data in HCI is complicated and requires significant pre-processing and post-processing to enable the discovery and characterization of exoplanets. Angular differential imaging (ADI) is an observational technique for HCI that utilizes the Earth's rotation throughout a night of observing [@liu:2004; @marois_angular_2006]. Normally telescopes have optics to counter this rotation, but disabling these optics and taking images throughout the night will give a sequence of frames with rotatinng astrophysical signal. The systematic signal, though, is from the telescope optics, so it does not appear to rotate like the sky, but will slowly vary over time. The sequence of frames are pre-processed then co-aligned and concatenated together into a data cube. ADI algorithms exploit the difference in rotation to approximate and subtract the systematics. The remaining sequence of residual images can then be reoriented and combined to average out the noise and reveal potential companions.

The difference in ADI algorithms is primarily in how the systematic signal is estimated. These algorithms can be applied in different ways, though- for example instead of estimating the systematics from the *target* data cube, a *reference* can be used (reference differential imaging, RDI). Instead of using the entire frame, the cube can be processed in annuli corresponding to different rates of apperant rotation. These applications are separate from the underlying description of the algorithms. This is the same when considering GPU programming or out-of-core processing, which appear like implementation details in comparsion to the algorithm. Creating a *modular* and *generic* framework for HCI enables much more thorough and deeper investigations into the capabilities of HCI in the search for exoplanets.

# Statement of need

`ADI.jl` is a Julia framework for post-processing high-contrast imaging (HCI) data. By organizing algorithms separately from their application, `ADI.jl` offers a modular API the benefits both observers and algorithm designers. Observers can rapidly experiment with access to as much or as little of the underlying API to optimally reduce their data. Julia's dynamic just-in-time (JIT) LLVM compiler [@Julia-2017; @2018arXiv180803370B] means this experimentation comes at a low runtime cost to the observer, enabling broader experimentation or higher throughput, for example, in the case of survey pipelines.

Algorithm designers will find that Julia is highly composable, so extending or adding a new algorithm only requires writing the code that is *unique to that algorithm*. Julia's language interoperability also means the algorithm can be implemented in Python or C, for example. The low barrier for entry combined with the growing feature set of `ADI.jl` makes writing new algorithms easy and rewarding.

Currently `ADI.jl` supports full-frame ADI and RDI reduction, with experimental support for SDI reduction. The algorithms that are currently implemented are median subtraction [@marois_angular_2006], principal component analysis [PCA/KLIP; @soummer_detection_2012], non-negative matrix factorization [NMF; @ren_non-negative_2018], and fixed-point greedy disk subtraction [GreeDS; @pairet_reference-less_2019; @pairet_mayonnaise_2020]. In addition, common metrics such as S/N maps and contrast curves are available for posterior analysis. Forward modeling is being built in a separate Julia package, `Firefly.jl`, as part of active research.


# Comparisons with existing software

High-contrast imaging as a field predominantly utilizes Python for data reduction. We break down some of the necessary computations into *pre-processing*, which includes raw calibration of data, the centering and stacking of the data cube, bad-pixel removal, etc., *post-processing*, which includes the PSF approximation and subtraction, *posterior analysis* which includes analyzing the results of post-processing like detection statistics and throughput calibration, and finally *forward modeling* which includes various statistical models for companions and disks for use with post-processing algorithms. `ADI.jl` primarily focuses on post-processing and posterior analysis.

Some notable libraries for HCI tasks include the Vortex Imaging Pipeline (`VIP`) [@gomez_gonzalez_vip_2016], `pyKLIP` [@2015ascl.soft06001W], and `PynPoint` [@pynpoint:2019]. A table of the feature sets of these packages alongside `ADI.jl` is presented below. In particular, `VIP` has served as a useful source of information regarding HCI image-processing as well as detailed implementations of common ADI algorithms. This has been indispensable in the development of `ADI.jl`, although this package is not a direct translation.

In general `VIP` offers the most diversity in algorithms and their applications, but not all algorithms are as feature-complete as the PCA implementation. `VIP` also contains many useful utilizies for pre-processing and a pipeline framework. `pyKLIP` only uses the PCA algorithm, but offers many more forward modeling implementations than any other package. `PynPoint` has a highly modular pre-processing module that is focused on pipelines.


Table 1: Comparison of features across different HCI frameworks. Techniques marked with * indicate partial support, meaning that not all algorithms are supported

 Framework | Pre-processing | Algorithms | Techniques | Posterior analysis | Forward modeling
-|-|-|-|-|-
ADI.jl | no | median, PCA, NMF, fixed-point GreeDS | Full-frame ADI/RDI, SDI (experimental) | detection maps, contrast curve | no
VIP | yes | median, LOCI, PCA, NMF, LLSG, ANDROMEDA | Full-frame ADI/RDI, SDI, annular ADI/RDI* | detection maps, blob detection, STIM, ROC, contrast curve | NegFC
pyKLIP | no | PCA | Full-frame ADI/RDI, SDI, annular ADI/RDI | detection maps, blob detection, contrast curve | KLIP-FM, Planet Evidence, matched filter (FMMF), spectrum fitting, DiskFM
PynPoint | yes | median, PCA | Full-frame ADI/RDI, SDI | detection maps, contrast curve | no

# Acknowledgments

We acknowledge the `VIP` team for creating a package that has been an indispensable resource in the creation of `ADI.jl`. We also acknowledge the ongoing research using `ADI.jl` for forward modeling of exoplanets with the package `Firefly.jl`.

# References
