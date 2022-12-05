#' YlGnBu palette.
#'
#' YlGnBu
#' @description Create a YlGnBu palette.
#' @param paletteLength.num <numeric>: color number.
#' @param space.chr <numeric>: a character string; interpolation in RGB or CIE Lab color spaces. See ?grDevices::colorRamp for more details. (Default "rgb")
#' @param interpolate.chr <numeric>: use spline or linear interpolation. See ?grDevices::colorRamp for more details. (Default "linear")
#' @param bias.num <numeric>: a positive number.  Higher values give more widely spaced colors at the high end. See ?grDevices::colorRamp for more details. (Default 1)
#' @return  A vector of color.
#' @examples
#' YlGnBu(9)
#'
YlGnBu <- function(paletteLength.num = NULL, space.chr = "rgb", interpolate.chr = "linear",
                   bias.num = 1) {
    (grDevices::colorRampPalette(
        colors = c(
            "#FFFFD9", "#EDF8B1", "#C7E9B4",
            "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58"
        ),
        space = space.chr, interpolate = interpolate.chr, bias = bias.num
    ))(paletteLength.num)
}
