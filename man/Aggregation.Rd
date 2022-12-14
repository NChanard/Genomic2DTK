% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Aggregation.R
\name{Aggregation}
\alias{Aggregation}
\title{Aggregation of matrices list.}
\usage{
Aggregation(
  ctrlMatrices.lst = NULL,
  matrices.lst = NULL,
  minDist.num = NULL,
  maxDist.num = NULL,
  trans.fun = NULL,
  agg.fun = "mean",
  rm0.bln = FALSE,
  diff.fun = "substraction",
  scaleCorrection.bln = FALSE,
  correctionArea.lst = NULL,
  statCompare.bln = FALSE,
  orientate = TRUE
)
}
\arguments{
\item{ctrlMatrices.lst}{<list[matrix]>: The matrices list to aggregate as control.}

\item{matrices.lst}{<list[matrix]>: The matrices list to aggregate.}

\item{minDist.num}{<numeric>: The minimal distance between anchor and bait.}

\item{maxDist.num}{<numeric>: The maximal distance between anchor and bait.}

\item{trans.fun}{<function or chracter>: The function use to transforme or scale values in each submatrix before aggregation. If the parameter is character so:
\itemize{
\item "quantile" or "qtl" apply function dplyr::ntile(x,500)
\item "percentile" or "prct" apply percentile.
\item "rank" apply a ranking.
\item "zscore" apply a scaling.
\item "minmax" apply a Genomic2DTK::MinMaxScale.
\item "mu" apply a Genomic2DTK::MeanScale.
\item other or NULL don't apply transformation (Default).
}}

\item{agg.fun}{<function or chracter>: The function use to aggregate each pixel in matrix. If the parameter is a character so:
\itemize{
\item "50%" or "median" apply the median
\item "+" or "sum" apply the sum
\item other (Default) apply the mean
}}

\item{rm0.bln}{<logical>: Whether 0 are replace with NA. (Default FALSE)}

\item{diff.fun}{<function or chracter>: The function use to compute differential. If the parameter is character so:
\itemize{
\item "-", "substract" or "substraction" apply a substraction
\item "/" or "ratio" apply a ratio
\item "log2","log2-","log2/" or "log2ratio" apply a log2 on ratio
\item other (Default) apply a log2 on 1+ratio
}}

\item{scaleCorrection.bln}{<logical>: Whether a correction should be done on the median value take in ane noising area. (Default TRUE)}

\item{correctionArea.lst}{<list>: Nested list of indice that define a noising area fore correction. List muste contain in first an element "i" (row indices) then an element called "j" (columns indices). If NULL automatically take in upper left part of aggregated matrices. (Default NULL)}

\item{statCompare.bln}{<logical>: Whether a t.test must be apply to each pxl of the differential aggregated matrix.}

\item{orientate}{<logical>: Whether matrices must be orientate before the aggregation.}
}
\value{
A matrix
}
\description{
Aggregates all the matrices of a list (or two lists in case of differential aggregation) into a single matrix. This function allows to apply different aggregation (average, sum, ...), transformation (rank, percentage, ...) and differential (subtraction, ratio, ...) functions.
}
\details{
Aggregation
}
\examples{
# Data
data(Beaf32_Peaks.gnr)
data(HiC_Ctrl.cmx_lst)
data(HiC_HS.cmx_lst)

# Index Beaf32
Beaf32_Index.gnr <- IndexFeatures(
    gRange.gnr_lst = list(Beaf = Beaf32_Peaks.gnr),
    chromSize.dtf = data.frame(seqnames = c("2L", "2R"), seqlengths = c(23513712, 25286936)),
    binSize.num = 100000
)

# Beaf32 <-> Beaf32 Pairing
Beaf_Beaf.gni <- SearchPairs(indexAnchor.gnr = Beaf32_Index.gnr)
Beaf_Beaf.gni <- Beaf_Beaf.gni[seq_len(2000)] # subset 2000 first for exemple

# Matrices extractions center on Beaf32 <-> Beaf32 point interaction
interactions_Ctrl.mtx_lst <- ExtractSubmatrix(
    feature.gn = Beaf_Beaf.gni,
    hic.cmx_lst = HiC_Ctrl.cmx_lst,
    referencePoint.chr = "pf"
)
interactions_HS.mtx_lst <- ExtractSubmatrix(
    feature.gn = Beaf_Beaf.gni,
    hic.cmx_lst = HiC_HS.cmx_lst,
    referencePoint.chr = "pf"
)

# Aggregate matrices in one matrix
aggreg.mtx <- Aggregation(interactions_Ctrl.mtx_lst)

# Differential Aggregation
aggregDiff.mtx <- Aggregation(
    ctrlMatrices.lst = interactions_Ctrl.mtx_lst,
    matrices.lst = interactions_HS.mtx_lst
)

}
