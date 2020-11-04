---
title: 'ADI.jl: A Julia Suite for High-Contrast Imaging'
tags:
  - Julia
  - astronomy
  - high-contrast imaging
  - direct imaging
  - image processing
authors:
  - name: Miles D. Lucas
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

# Acknowledgements

# References
