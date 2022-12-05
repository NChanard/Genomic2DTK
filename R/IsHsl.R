#' Check HSL color format.
#'
#' IsHsl
#' @description Check if a color is in HSL color format.
#' @param color.col <character or numeric>: a color.
#' @return a logical.
#' @examples
#' IsHsl("red")
#' IsHsl("#FFFFFF")
#' IsHsl(c(125, 125, 125))
#' IsHsl(c(43.8, 0.873, 0.492))
#'
IsHsl <- function(
    color.col = NULL
) {
    logical.bln <- lapply(
            color.col[-1],
            function(value.num) {0 <= value.num & value.num <= 1 }
        ) |>
        unlist() |>
        sort()
    logical.bln <- logical.bln[[1]]
    return(
        class(color.col) %in% c("list", "numeric", "integer") &&
        0 <= color.col[[1]] &&
        color.col[[1]] < 360 &&
        logical.bln
    )
}
