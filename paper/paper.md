---
title: 'ADI.jl: A Julia Suite for High-Contrast Imaging'
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

High-contrast imaging (HCI) is a powerful technique for discovering and characterizing exoplanets. Being able to probe the architecture, formation, and atmospheres of planets *directly* is necessary for advancing companion formation and evolution theory. The process required to image a planet is daunting, however, due to the brightness of proximity of their host stars. One part of making such a difficult detection is the image processing of HCI data to remove systematic signal and attenuate noise. The size of HCI data and complexity of processing algorithms requires an efficient numerical framework while simultaneously offering modularity for rapid experimentation and discovery.

# Statement of need

`ADI.jl` is a Julia framework for the post-processing of HCI data. The high-level composability and the low-level efficiency of Julia is uniquely suited for HCI image-processing. By organizing the process of the data reduction into a common interface in `ADI.jl`, implementing a new algorithm requires minimal code while offering a maximal API. This dynamic combination is designed with both observationalists and algorithm creators in mind. By having an efficient and modular system observationalists can quickly try different algorithms and application methods for the rapid discovery and characterization of exoplanets. Similarly, the composable API is naturally extensible and new algorithms or extensions of existing algorithms require a minimal implementation while providing the full API that `ADI.jl` has to offer.

# A common framework for HCI

The data in high-contrast imaging is complicated and requires significant pre-processing and post-processing to enable the discovery and characterization of exoplanets. In particular, the post-processing requires advanced algorithms for removing systematic signal and attenuating noise. Fortunately, though, many of these algorithms follow a similar process for reducing HCI data.

For a brief primer, angular differential imaging (ADI) is an observing technique that utilizes the Earth's rotation over a night of observing. The data it delivers are a sequence of images over time, often stacked into a data cube. This cube is comprised of systematics and astrophysical signal, but due to the Earth's rotation the astrophysical signal appears to rotate throughout the sequence. ADI algorithms exploit this difference in structure to approximate and subtract the systematics. The remaining sequence of residual images can then be reoriented and combined to attenuate the noise and reveal potential companions.

In `ADI.jl` this process is implemented in a functional framework. Each algorithm implements a `reconstruct` method, which describes the approximation of the systematic signal given a *target* data cube, a vector of parallactic angles (if needed), and can optionally use a separate *reference* data cube. The algorithm is a first-class function, which means the options specific to the algorithm can be separated from the arguments for the usage of the algorithm. For example, the common reduction algorithm of principla components analysis (PCA) requires specifying the rank of the low-rank basis. This parameter is set within the `PCA` algorithm which is then used generically with the `ADI.jl` frameowrk.

The `reconstruct` method is all that is needed to describe an algorithm since the reorientation and collapsing of the target data cube is a common pattern. As such, applying the algorithm as a function will combine the reconstruction with the remainder of the reduction. This API is enabled by Julia's multiple-dispatch, which means functions dispatch on the type of all input arguments. This is in contrast to common single-dispatch languages, like C and Python, which can only organize their algorithms using inheritance or function overloading^[function overloading is distinct from multiple-dispatch because function overloading only dispatches on the overloaded caller, not on any of the other types. This is roughly the same reason why mix-ins or monkeypatching are not multiple-dispatch.]. As a result, the Julia ecosystem is much more composable, since there is no conflict with extending a function for a new type (or conversely, a new function to extend a type).

The composability that multiple-dispatch offers is beneficial to observers and algorithm designers. Observers will find that treating algorithms as building blocks allows experimenting with different reductions is quicker and easier, which leads to more rapid exploration of hyperparameters leading to more powerful and accurate science. Julia's dynamic nature and ecosytem of tools designed for discovery (such as IJulia, weave, Pluto, and VSCode) is also well suited for observers. Algorithm designers exploit composability by only having write the core essence of their algorithm and having it fit in to the existing framework seamlessly. In addition, for algorithms that do not fit into the precise building blocks can always override the default implementations. Lastly, Julia's rich language inter-operability allows using algorithms defined in Python, C, Fortran, or R, for example, creating an HCI framework that can be useful to researchers despite different programming needs.


# Comparison with existing software

High-contrast imaging, as a field, predominantly utilizes Python for data reduction. We break down some of the necesssary computaions into *pre-processing*, which includes raw calibration of data, the centering and stacking of the data cube, bad-pixel removal, etc., *post-processing*, which includes the PSF approximation and subtraction, *posterior analysis* which includes analyzing the results of post-processing like detection statistics and throughput calibration, and finally *forward modeling* which includes various statistical models for companions and disks for use with post-processing algorithms.

Some notable libraries for HCI tasks include the Vortex Imaging Pipeline (`VIP`) **citation**, `pyKLIP` **citation**, and `PynPoint` **citation**. A table of the feature sets of these packages alongside `ADI.jl` is presented in **table link**. In particular, VIP has served as a useful source of information regarding HCI image-processing as well as detailed implementations of common ADI algorithms. This has been indispensible in the development of `ADI.jl`, although this package is not a direct translation.

VIP offers the most diversity in algorithms and their applications, but at the cost of flexibility. The implementations of the algorithms suffer from long lists of arguments, large blocks of parsing the input arguments and coercing them into usable types or manually dispatching. Because of this, algorithms with similar techniques, such as annular processing using the median vs. annular processing with PCA, are not guaranteed to be applied the same way and any new algorithms must repeat code and add complexity in order to use the annular processing technique.

 Framework | Pre-processing | Algorithms | Techniques | Posterior analysis | Forward modeling 
-|-|-|-|-|-
ADI.jl | no | median, PCA, NMF, fixed-point GreeDS | Full-frame ADI/RDI, SDI (experimental) | detection maps, contrast curve | no ^[separate package in development, `Firefly.jl`]
VIP | yes | median, LOCI, PCA, NMF, LLSG, ANDROMEDA | Full-frame ADI/RDI, SDI, annular ADI/RDI for some algorithms | detection maps, blob detection, STIM, ROC, contrast curve | NegFC 
pyKLIP | no | PCA | Full-frame ADI/RDI, SDI, annular ADI/RDI | detection maps, blob detection, contrast curve | KLIP-FM, Planet Evidence, matched filter (FMMF), spectrum fitting, DiskFM
PynPoint | yes | median, PCA | Full-frame ADI/RDI, SDI | detection maps, contrast curve | no 


# Acknowledgements

# References
