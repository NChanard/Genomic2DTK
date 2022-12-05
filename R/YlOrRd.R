#' YlOrRd palette.
#'
#' YlOrRd
#' @description Create a YlOrRd palette.
#' @param paletteLength.num <numeric>: color number.
#' @param space.chr <numeric>: a character string; interpolation in RGB or CIE Lab color spaces. See ?grDevices::colorRamp for more details. (Default "rgb")
#' @param interpolate.chr <numeric>: use spline or linear interpolation. See ?grDevices::colorRamp for more details. (Default "linear")
#' @param bias.num <numeric>: a positive number. Higher values give more widely spaced colors at the high end. See ?grDevices::colorRamp for more details. (Default 1)
#' @return A vector of color.
#' @examples
#' YlOrRd(9)
#'
YlOrRd <- function(
    paletteLength.num = NULL, space.chr = "rgb",
    interpolate.chr = "linear", bias.num = 1
) {
    (grDevices::colorRampPalette(
        colors = c(
            "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C",
            "#FC4E2A", "#E31A1C", "#BD0026", "#800026"
        ),
        space = space.chr,
        interpolate = interpolate.chr,
        bias = bias.num
    ))(paletteLength.num)
}
