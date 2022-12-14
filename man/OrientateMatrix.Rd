% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OrientateMatrix.R
\name{OrientateMatrix}
\alias{OrientateMatrix}
\title{Matrix orientation}
\usage{
OrientateMatrix(matrice.mtx)
}
\arguments{
\item{matrice.mtx}{<matrix or List[matrix]>: Matrix or matricies list to oriente}
}
\value{
Oriented matrix or matricies list
}
\description{
Oriente extracted Matrix according to the anchors and bait order. Apply a 180° rotation follow with a transposation on a matrix or on matricies in a list according to the interactions attributes of the list.
}
\details{
OrientateMatrix
}
\examples{
# Data
data(Beaf32_Peaks.gnr)
data(HiC_Ctrl.cmx_lst)

# Index Beaf32 in TADs domains
Beaf32_Index.gnr <- IndexFeatures(
    gRange.gnr_lst = list(Beaf = Beaf32_Peaks.gnr),
    chromSize.dtf = data.frame(
        seqnames = c("2L", "2R"),
        seqlengths = c(23513712, 25286936)
    ),
    binSize.num = 100000
)

# Beaf32 <-> Beaf32 Pairing
Beaf_Beaf.gni <- SearchPairs(indexAnchor.gnr = Beaf32_Index.gnr)
Beaf_Beaf.gni <- Beaf_Beaf.gni[seq_len(2000)] # subset 2000 first for exemple

# Matrices extractions center on Beaf32 <-> Beaf32 point interaction
interactions_PF.mtx_lst <- ExtractSubmatrix(
    feature.gn = Beaf_Beaf.gni,
    hic.cmx_lst = HiC_Ctrl.cmx_lst,
    referencePoint.chr = "pf"
)

# Matrices Orientation
oriented_Interactions_PF.mtx_lst <- OrientateMatrix(interactions_PF.mtx_lst)

}
