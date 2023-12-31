% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{getChannelSummaryStats}
\alias{getChannelSummaryStats}
\title{Calculate a summary statistic of some channels of
all flowFrames of a flowSet}
\usage{
getChannelSummaryStats(
  fs,
  channels = NULL,
  statFUNs = stats::median,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{fs}{a flowCore::flowSet}

\item{channels}{which channels (integer index(ices) or character(s)):
\itemize{
\item if it is a character vector,
it can refer to either the channel names, or the marker names
\item if it is a numeric vector,
it refers to the indexes of channels in \code{fs}
\item if NULL all scatter and fluorescent channels of \code{fs}
will be selected
}}

\item{statFUNs}{a list (possibly of length one) of
functions to call to calculate the statistics, or a simple function
This list can be named, in that case, these names will be transfered to the
returned value.}

\item{verbose}{if \code{TRUE}, output a message
after each single distance calculation}

\item{...}{additional parameters passed to \code{getEMDDist()}}
}
\value{
a list of named statistic matrices.
In each stat matrix, the columns are the channel statistics
for all flowFrames of the flowSet.
Exception: if only one stat function (and not a list) is passed in
\code{statFUNs}, the return value is simplified to the stat matrix itself.
}
\description{
Calculate a summary statistic of some channels of
all flowFrames of a flowSet
}
\examples{

library(CytoPipeline)

data(OMIP021Samples)

# estimate scale transformations 
# and transform the whole OMIP021Samples

transList <- estimateScaleTransforms(
    ff = OMIP021Samples[[1]],
    fluoMethod = "estimateLogicle",
    scatterMethod = "linearQuantile",
    scatterRefMarker = "BV785 - CD3")

OMIP021Trans <- CytoPipeline::applyScaleTransforms(
    OMIP021Samples, 
    transList)

channelsOrMarkers <- c("FSC-A", "SSC-A", "BV785 - CD3")

# calculate mean for each 4 selected channels, for each 2 samples

channelMeans <- getChannelSummaryStats(
    OMIP021Trans,
    channels = channelsOrMarkers,
    statFUNs = mean)
    
# calculate median AND std deviation
# for each 4 selected channels, for each 2 samples

channelMedians <- getChannelSummaryStats(
    OMIP021Trans,
    channels = channelsOrMarkers,
    statFUNs = list("median" = stats::median, 
                    "std.dev" = stats::sd))
   
}
