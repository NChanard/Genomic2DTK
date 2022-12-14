#' Find threshold for outliers based on quantiles.
#'
#' QtlThreshold
#' @keywords internal
#' @description Find threshold for outliers triming based on quantiles.
#' @param x.num <numeric>: Numeric vector.
#' @param prct.num <numeric>: Percentage (0-100) threshold. (Default 5)
#' @param bounds.chr <character>: Bounds to return, "lower", "upper" or "both". (Default "both")
#' @return Numerical vector of thresholds values for outliers triming.
#' @examples
#' set.seed(1111)
#' x.num <- 0:100
#' x.num <- sort(x.num)
#' x.num
#' QtlThreshold(x.num, prct.num = 5, bounds.chr = "lower")
#' QtlThreshold(x.num, prct.num = 5, bounds.chr = "both")
#' QtlThreshold(x.num, prct.num = 5, bounds.chr = "upper")
#'
QtlThreshold <- function(
    x.num = NULL, prct.num = 5, bounds.chr = "both"
) {
    probs.num <- dplyr::case_when(
        bounds.chr == "both" ~ c(prct.num / 200, 1 - (prct.num / 200)),
        bounds.chr == "upper" ~ c(NA, 1 - (prct.num / 100)),
        bounds.chr == "lower" ~ c(prct.num / 100, NA)
    )
    return(stats::quantile(x.num, na.rm = TRUE, probs.num))
}
