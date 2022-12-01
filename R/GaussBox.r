#' One dimension Gaussian kernel.
#'
#' GaussBox
#' @description One dimension Gaussian kernel.
#' @param sd.num <numeric>: standard deviation parameter of the gaussian. (Default 1)
#' @param boxSize.num <numeric>: size of kernel. If NULL size is 1+4*sd.num (Default NULL)
#' @param scale.chr <character>: scaling kind of box. If "1" sum of kernel equal 1. If "int" Minimal value of kernel is 1 and all entry are integer. If "none", kernel is not scale. (Default "1")
#' @return numerical vector.
#' @examples
#' GaussBox(sd.num=5, scale.chr="none")
#' GaussBox(scale.chr="1")
#' GaussBox(scale.chr="int")

GaussBox <- function(sd.num = 1, boxSize.num = NULL, scale.chr = "1") {
    if (is.null(boxSize.num)) {
        boxSize.num <- 1 + 4 * sd.num
    }
    x.num <- as.vector(scale(seq_len(boxSize.num), scale = FALSE, center = TRUE))
    box <- lapply(x.num, function(x) {
        xInterval.num <- seq((x - 0.5), (x + 0.5), by = 0.01) |>
            lapply(function(xi) {
                Gauss(x = xi, sd.num = sd.num)
            }) |>
            unlist() |>
            mean()
        return(xInterval.num)
    }) |>
        unlist()
    if (scale.chr == "1") {
        box <- {
            box/sum(abs(box))
        }
    } else if (scale.chr == "int") {
        box <- ceiling(box/box[1])
    }
    return(box)
}