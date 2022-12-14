% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICEnorm.R
\name{ICEnorm}
\alias{ICEnorm}
\title{Compute Iterative Correction.}
\usage{
ICEnorm(hic.cmx, qtlTh.num = 0.15, maxIter.num = 50)
}
\arguments{
\item{hic.cmx}{<contactMatrix>: The HiC maps chunk to normalize.}

\item{qtlTh.num}{<numerical>: The threshold quantile below which the bins will be ignored. (Default 0.15)}

\item{maxIter.num}{<numerical>: The maximum iteration number.}
}
\value{
A normalized contactMatrix
}
\description{
Compute Iterative Correction (Vanilla Count) on hic maps.
}
\details{
ICEnorm
}
\examples{
data(HiC_Ctrl.cmx_lst)
HiC_Ctrl_ICE.cmx <- ICEnorm(HiC_Ctrl.cmx_lst[['2L_2L']])

}
