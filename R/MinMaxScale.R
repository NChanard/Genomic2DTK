#' Scales values on min-max range.
#'
#' MinMaxScale
#' @description Scale values on min-max range.
#' @param x.num <numeric>: numerical vector.
#' @param min.num <numeric>: minimal value after scaling.
#' @param max.num <numeric>: maximal value after scaling.
#' @return scaled numeric vector.
#' @examples
#' set.seed(655213)
#' x.num <- rnorm(500, 500)
#' set.seed(522613)
#' y.num <- rnorm(500, 100)
#' plot(density(x.num), col = "red", xlim = c(min(y.num), max(x.num)))
#' lines(density(y.num), col = "green")
#' plot(density(MinMaxScale(x.num)), col = "red", xlim = c(min(MinMaxScale(y.num)), max(MinMaxScale(x.num))))
#' lines(density(MinMaxScale(y.num)), col = "green")
#'
MinMaxScale <- function(x.num, min.num = (0), max.num = 1) {
    min.num + (((x.num - min(x.num, na.rm = TRUE)) * (max.num - min.num)) / (max(x.num,
        na.rm = TRUE
    ) - min(x.num, na.rm = TRUE)))
}
