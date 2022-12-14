% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rgb2Hex.R
\name{Rgb2Hex}
\alias{Rgb2Hex}
\title{Convert RGB to Hex.}
\usage{
Rgb2Hex(rgb.col = NULL, alpha.bln = FALSE)
}
\arguments{
\item{rgb.col}{<integer>: An integer of the color's RGB code.}

\item{alpha.bln}{<logical>: Whether the alpha layer should be returned. (Default FALSE)}
}
\value{
A character of the color's hexadecimal code.
}
\description{
Convert a color in RGB format to hexadecimal format.
}
\details{
Rgb2Hex
}
\examples{
Rgb2Hex(c(235, 176, 16, 127), alpha.bln = TRUE)

}
\keyword{internal}
