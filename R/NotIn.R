#' Exclusion binary operator.
#'
#' NotIn
#' @keywords internal
#' @description Binary operator, inverse to \%in\%.
#' @param lhs <vector or NULL>: Values to be compared against rhs
#' @param rhs <vector or NULL>: Values to be compared against lhs
#' @return A boolean.
#' @examples
#' "A" |> NotIn(c("A", "B", "C"))
#' "A" |> NotIn(c("B", "C", "D"))
#' NotIn("A", c("A", "B", "C"))
#' NotIn("A", c("B", "C", "D"))
#'
"NotIn" <- function(
    lhs, rhs
) {
    return(match(lhs, rhs, nomatch = 0L) == 0L)
}
