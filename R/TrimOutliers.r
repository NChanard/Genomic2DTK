#' Remove outliers.
#'
#' TrimOutliers
#' @description Replace values of a numerical vector that are below a minimal thresholds and/or above maximal thresholds.
#' @param x.num <numeric>: numeric vector.
#' @param tresholds.num <numeric>: numeric vector of length 2. first value is minimal threshold, second value maximal threshold (Default find threshold based on standarrd deviation. see `SdThreshold` function)
#' @param clip.bln <logical>: if TRUE the value out of bounds are replace with threshodls values. If FALSE the Values out of bound are replace with NA (Default FALSE).
#' @return trimed Numerical vector.
#' @examples
#' set.seed(1111)
#' x.num = rnorm(1000)
#' x.num = sort(x.num)
#' x.num[990:1000]
#' SdThreshold(x.num)
#' TrimOutliers(x.num)[990:1000]
#' TrimOutliers(x.num, clip=TRUE)[990:1000]

TrimOutliers <- function(x.num, tresholds.num=SdThreshold(x.num), clip.bln=FALSE){
        if( clip.bln){
            x.num[which(x.num>tresholds.num[2])] <-  tresholds.num[2]
            x.num[which(x.num<tresholds.num[1])] <-  tresholds.num[1]
        }else{
            x.num[which(x.num>tresholds.num[2] | tresholds.num[1]>x.num) ] <- NA
        }
        return(x.num)
}