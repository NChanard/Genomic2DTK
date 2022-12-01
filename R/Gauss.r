#' Gaussian formula.
#'
#' Gauss
#' @description Gaussian formula in 1 or 2 dimension.
#' @param x <numeric>: x value.
#' @param y <numeric>: y value for 2 dimensional gaussian.
#' @param sd.num <numeric>: standard deviation parameter of the gaussian. (Default 1)
#' @param mu <numeric>:  mean deviation parameter of the gaussian. (Default 0)
#' @return Result of Gaussian formula
#' @examples
#' Gauss(x=1)
#' Gauss(x=1,y=2)

Gauss <- function(x = NULL, y = NULL, sd.num = 1, mu = 0) {
    x <- x[1]
    y <- y[1]
    if (is.null(y)) {
        return(1/(sd.num * sqrt(2 * pi)) * exp(-((x - mu)^2)/(2 * sd.num^2)))
    } else {
        return(1/(2 * pi * sd.num^2) * exp(-((x^2 + y^2)/(2 * sd.num^2))))
    }
}