#' Convert RGB to Hex.
#'
#' Rgb2Hex
#' @keywords internal
#' @description Convert a color in RGB format to hexadecimal format.
#' @param rgb.col <integer>: An integer of the color's RGB code.
#' @param alpha.bln <logical>: Whether the alpha layer should be returned. (Default FALSE)
#' @return A character of the color's hexadecimal code.
#' @examples
#' Rgb2Hex(c(235, 176, 16, 127), alpha.bln = TRUE)
#'
Rgb2Hex <- function(
    rgb.col = NULL, alpha.bln = FALSE
) {
    if (3 > length(rgb.col) | length(rgb.col) > 4) {
        err.chr <- paste0(
            "Need 3 or 4 values beetween 0 and 255, first value for red, ",
            "second for green, third for blue and last for alpha"
        )
        stop(err.chr)
    } else {
        if (IsRgb(rgb.col)) {
            if (length(rgb.col) == 3) {
                rgb.col <- c(rgb.col, 255)
            }
            hex.col <- lapply(rgb.col, function(value.num) {
                fisrtBit <- value.num %/% 16
                if (fisrtBit > 9) {
                    fisrtBit <- letters[fisrtBit - 9]
                }
                secondBit <- value.num %% 16
                if (secondBit > 9) {
                    secondBit <- letters[secondBit - 9]
                }
                return(c(fisrtBit, secondBit))
            }) |>
                unlist()
            hex.col <- paste0(c("#", hex.col), collapse = "")
        } else {
            err.chr <- paste0(
                "Need 3 or 4 values beetween 0 and 255, first value for red, ",
                "second for green, third for blue and last for alpha"
            )
            stop(err.chr)
        }
    }
    if (!alpha.bln) {
        hex.col <- substr(hex.col, 1, nchar(hex.col) - 2)
    }
    return(hex.col)
}
