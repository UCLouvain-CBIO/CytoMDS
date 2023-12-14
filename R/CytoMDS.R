# CytoMDS - Copyright (C) <2023>
# <Université catholique de Louvain (UCLouvain), Belgique>
#
#   Description and complete License: see LICENSE file.
#
# This program (CytoMDS) is free software:
#   you can redistribute it and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details (<http://www.gnu.org/licenses/>).


#' @title CytoMDS package
#' 
#' @docType package
#' 
#' @name CytoMDS
#'
#' @rdname CytoMDS
#' 
#' @seealso [ggplotSampleMDS], [ggplotSampleMDSShepard], [computeMetricMDS]
#'
#' @description
#'
#' `CytoMDS` implements a low dimensional visualization of a set of cytometry 
#' samples, in order to visually assess the 'distances' between them.
#' This, in turn, can greatly help the user to identify quality issues 
#' like batch effects or outlier samples, and/or check the presence of 
#' potential sample clusters that might align with the experimental design.  
#'   
#' The CytoMDS algorithm combines, on the one hand, 
#' the concept of Earth Mover's Distance (EMD), a.k.a. Wasserstein metric 
#' and, on the other hand, the Multi Dimensional Scaling (MDS) algorithm 
#' for the low dimensional projection.  
#' 
#' Also, the package provides some diagnostic tools for 
#' both checking the quality of the MDS projection, 
#' as well as tools to help with the interpretation of 
#' the axes of the projection.
#' 
NULL