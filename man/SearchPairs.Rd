% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SearchPairs.r
\name{SearchPairs}
\alias{SearchPairs}
\title{SearchPairs}
\usage{
SearchPairs(
  indexAnchor.gnr = NULL,
  indexBait.gnr = NULL,
  minDist.num = NULL,
  maxDist.num = NULL,
  verbose.bln = TRUE,
  cores.num = 1
)
}
\arguments{
\item{indexAnchor.gnr}{<GRanges>: GRanges used in the first part of the pairs.}

\item{indexBait.gnr}{<GRanges>: GRanges used in the second part of the pars. If NULL is replace with indexAnchor.gnr (Default NULL)}

\item{minDist.num}{<numeric>: the minimal distance between first and second parts of pairs. (Default NULL)}

\item{maxDist.num}{<numeric>: the maximal distance between first and second parts of pairs. (Default NULL)}

\item{verbose.bln}{<logical>: A logical value. If TRUE show the progression in console. (Default TRUE)}

\item{cores.num}{<integer> : An integer to specify the number of cores. (Default 1)}
}
\value{
A GInteractions object.
}
\description{
Compute pairs of genomic features that share a same range constraintes.
}
\examples{
library(GenomicED)
data(anchors_Index.gnr)
anchors_Index.gnr[1]
data(baits_Index.gnr)
baits_Index.gnr[1]
interactions.gni <- SearchPairs(
    indexAnchor.gnr = anchors_Index.gnr,
    indexBait.gnr   = baits_Index.gnr,
    minDist.num     = 9000, 
    maxDist.num     = 11000,
    cores.num       = 1,
    verbose.bln     = FALSE
    )
interactions.gni[1]
}