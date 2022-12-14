% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QtlThreshold.R
\name{QtlThreshold}
\alias{QtlThreshold}
\title{Find threshold for outliers based on quantiles.}
\usage{
QtlThreshold(x.num = NULL, prct.num = 5, bounds.chr = "both")
}
\arguments{
\item{x.num}{<numeric>: Numeric vector.}

\item{prct.num}{<numeric>: Percentage (0-100) threshold. (Default 5)}

\item{bounds.chr}{<character>: Bounds to return, "lower", "upper" or "both". (Default "both")}
}
\value{
Numerical vector of thresholds values for outliers triming.
}
\description{
Find threshold for outliers triming based on quantiles.
}
\details{
QtlThreshold
}
\examples{
set.seed(1111)
x.num <- 0:100
x.num <- sort(x.num)
x.num
QtlThreshold(x.num, prct.num = 5, bounds.chr = "lower")
QtlThreshold(x.num, prct.num = 5, bounds.chr = "both")
QtlThreshold(x.num, prct.num = 5, bounds.chr = "upper")

}
\keyword{internal}
