#' Check RGB color format.
#'
#' IsRgb
#' @description Check if a color is in RGB color format.
#' @param color.col <character or numeric>: a color.
#' @return a logical.
#' @examples
#' IsRgb("red")
#' IsRgb("#FFFFFF")
#' IsRgb(c(125, 125, 125))
#' IsRgb(c(43.8, 0.873, 0.492))
#'
IsRgb <- function(
    color.col = NULL
) {
    logical.bln <- lapply(
            color.col,
            function(value.num) {0 <= value.num & value.num <= 255}
        ) |>
        unlist() |>
        sort()
    logical.bln <- logical.bln[[1]]
    return(
        (!IsHsl(color.col)) &&
        (class(color.col) %in% c("list", "numeric", "integer") &&
        logical.bln)
    )
}
