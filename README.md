# Low Dimensional Projection of Cytometry Samples

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/UCLouvain-CBIO/CytoMDS/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/UCLouvain-CBIO/CytoMDS/actions?query=workflow%3AR-CMD-check-bioc)
[![license](https://img.shields.io/badge/license-GPL3.0-blue)](https://opensource.org/licenses/GPL-3.0)
[![codecov.io](https://codecov.io/github/UCLouvain-CBIO/CytoMDS/coverage.svg?branch=main)](https://codecov.io/github/UCLouvain-CBIO/CytoMDS?branch=main)

The `CytoMDS` package  implements a low dimensional visualization of a set
of cytometry samples, in order to visually assess the 'distances' between them.
This, in turn, can greatly help the user to identify quality issues 
like batch effects or outlier samples, and/or check the presence of potential 
sample clusters that might align with the experimental design.  

The `CytoMDS` algorithm combines, on the one hand, the concept of Earth Mover's 
Distance (EMD), a.k.a. Wasserstein metric and, on the other hand, 
the Multi Dimensional Scaling (MDS) algorithm for the low dimensional
projection.  

Also, the package provides some diagnostic tools for both checking the quality 
of the MDS projection, as well as tools to help with the interpretation of 
the axes of the projection.
