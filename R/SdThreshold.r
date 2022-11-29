#' Find threshold for outliers triming based on standard deviation.
#'
#' SdThreshold
#' @description Find threshold for outliers triming based on standard deviation.
#' @param x.num <numeric>: numeric vector.
#' @param sdThreshold.num <numeric>: number of standard deviation. (Default 3)
#' @param bounds.chr <character>: bounds to return, "lower", "upper" or "both". (Default "both")
#' @return numerical vector of thresholds values for outliers triming
#' @examples
#' set.seed(1111)
#' x.num <- rnorm(1000)
#' x.num <- sort(x.num)
#' x.num
#' SdThreshold(x.num, sdThreshold.num=2, bounds.chr="lower")
#' SdThreshold(x.num, sdThreshold.num=2, bounds.chr="both")
#' SdThreshold(x.num, sdThreshold.num=2, bounds.chr="upper")

SdThreshold <- function(x.num=NULL, sdThreshold.num=3, bounds.chr="both"){
    mu.num <- mean(x.num,na.rm=TRUE)
    sdev.num <- stats::sd(x.num,na.rm=TRUE)
    upper.num <- mu.num + (sdThreshold.num*sdev.num)
    lower.num <- mu.num - (sdThreshold.num*sdev.num)
    tresholds.num <- dplyr::case_when(
        bounds.chr =="both" ~ c(lower.num,upper.num),
        bounds.chr =="upper" ~ c(NA,upper.num),
        bounds.chr =="lower" ~ c(lower.num,NA)
    )
    return(tresholds.num)
}
