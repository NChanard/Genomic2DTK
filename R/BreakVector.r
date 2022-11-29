#' Cut a vector.
#'
#' BreakVector
#' @description Compute the n+1 breaks of a vector in a linear or density based way with the possibility to fix minimal, center and maximal values.
#' @param x.num <numeric>: numerical vector.
#' @param min.num <numeric>: minimal fixed value.
#' @param center.num <numeric>: center fixed value.
#' @param max.num <numeric>: maximal fixed value.
#' @param n.num <numeric>: number of tile (return n.num+1 breaks).
#' @param method.chr <character>: kind of breaking. "linear" or "density". (Default "linear")
#' @return Numerical vector of breaks.
#' @examples
#' set.seed(31415)
#' x.num <- rnorm(100,50,200)
#' BreakVector(x.num=x.num, n.num=9)

BreakVector <- function(x.num=NULL, min.num=NULL, center.num=NULL, max.num=NULL, n.num=10, method.chr="linear") {
    n.num <- n.num+1
    if(method.chr=="linear"){
        x.num <- x.num[which(!is.na(x.num))]
        if(is.null(min.num)){min.num <- min(x.num, na.rm=TRUE)}
        if(is.null(max.num)){max.num <- max(x.num, na.rm=TRUE)}
        if(is.null(center.num)){
            breaks.num <- seq(min.num, max.num,length.out = n.num)
        }else if(min.num<center.num & center.num<max.num ){
            breaks.num <- c(seq(min.num, center.num,length.out = n.num%/%2+1), seq(center.num,max.num,length.out = n.num%/%2+1))
        }else{
            center.num <- stats::median(x.num, na.rm=TRUE)
            breaks.num <- c(seq(min.num, center.num,length.out = n.num%/%2+1), seq(center.num,max.num,length.out = n.num%/%2+1))
        }
    }else if(method.chr=="density"){
        breaks.num <- stats::quantile(x.num, na.rm=TRUE, probs=seq(0, 1, length.out=n.num))
    }else{
        stop("Error, method.chr muste be one of 'linear' or 'density'.\n")
    }
    return(breaks.num[!duplicated(breaks.num)])
}