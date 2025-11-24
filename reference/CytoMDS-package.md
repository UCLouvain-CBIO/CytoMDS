# CytoMDS: Low Dimensions projection of cytometry samples

This package implements a low dimensional visualization of a set of
cytometry samples, in order to visually assess the 'distances' between
them. This, in turn, can greatly help the user to identify quality
issues like batch effects or outlier samples, and/or check the presence
of potential sample clusters that might align with the exeprimental
design. The CytoMDS algorithm combines, on the one hand, the concept of
Earth Mover's Distance (EMD), a.k.a. Wasserstein metric and, on the
other hand, the Multi Dimensional Scaling (MDS) algorithm for the low
dimensional projection. Also, the package provides some diagnostic tools
for both checking the quality of the MDS projection, as well as tools to
help with the interpretation of the axes of the projection.

## See also

Useful links:

- <https://uclouvain-cbio.github.io/CytoMDS>

- Report bugs at <https://github.com/UCLouvain-CBIO/CytoMDS/issues>

## Author

**Maintainer**: Philippe Hauchamps <philippe.hauchamps@uclouvain.be>
([ORCID](https://orcid.org/0000-0003-2865-1852))

Authors:

- Laurent Gatto <laurent.gatto@uclouvain.be>
  ([ORCID](https://orcid.org/0000-0002-1520-2268))

Other contributors:

- Dan Lin <dan.8.lin@gsk.com> \[contributor\]
