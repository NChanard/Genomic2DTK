#' Sum by removing NA.
#'
#' Plus
#' @keywords internal
#' @description Perform sum by removing NA. If all values are NA return NA instead 0.
#' @param x <numerical>: a numerical vector
#' @return  The sum of numbers.
#' @examples
#' Plus(c(1,2,3))
#' sum(c(1,2,3),na.rm=TRUE)
#' sum(c(1,2,3))
#' Plus(c(1,2,NA))
#' sum(c(1,2,NA),na.rm=TRUE)
#' sum(c(1,2,NA))
#' Plus(c(NA,NA,NA))
#' sum(c(NA,NA,NA),na.rm=TRUE)
#' sum(c(NA,NA,NA))
Plus <- function(x){
    if(all(is.na(x))){
        c(x[0],NA)
    }else{
        sum(x,na.rm = TRUE)
    }
}