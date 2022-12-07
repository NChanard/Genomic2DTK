#' Add a value around a matrix.
#'
#' PadMtx
#' @description Add a value around a matrix.
#' @param mat.mtx <matrix>: Numerical matrix.
#' @param padSize.num <numeric>: Number of columns or rows to add. (Default 1)
#' @param value.num <numeric>: Value to add. If Null create mirror of choosen sides. (Default 0)
#' @param side.chr <character>: Side to pad, must be one or some of 'top','bot','right' or 'left'. (Default c('top','bot','right','left') )
#' @return A matrix.
#' @examples
#' mat.mtx <- matrix(seq_len(25), 5, 5)
#' PadMtx(
#'     mat.mtx = mat.mtx,
#'     padSize.num = 1,
#'     value.num = 0,
#'     side.chr = c("top", "bot", "right", "left")
#' )
#' PadMtx(
#'     mat.mtx = mat.mtx,
#'     padSize.num = 1,
#'     value.num = NULL,
#'     side.chr = c("top", "bot", "right", "left")
#' )
#' PadMtx(
#'     mat.mtx = mat.mtx,
#'     padSize.num = 1,
#'     value.num = 0,
#'     side.chr = c("right", "left")
#' )
#' PadMtx(
#'     mat.mtx = mat.mtx,
#'     padSize.num = 1,
#'     value.num = 0,
#'     side.chr = c("top")
#' )
#'
PadMtx <- function(
    mat.mtx = NULL, padSize.num = 1, value.num = 0,
    side.chr = c("top", "bot", "right", "left")
) {
    if ("top" %in% side.chr) {
        if (!is.null(value.num)) {
            row.lst <- rep(list(rep(value.num, dim(mat.mtx)[2])), padSize.num)
            row.pad <- do.call(rbind, row.lst)
        } else {
            row.pad <- mat.mtx[padSize.num:1, ]
        }
        mat.mtx <- rbind(row.pad, mat.mtx)
    }
    if ("bot" %in% side.chr) {
        if (!is.null(value.num)) {
            row.lst <- rep(list(rep(value.num, dim(mat.mtx)[2])), padSize.num)
            row.pad <- do.call(rbind, row.lst)
        } else {
            row.pad <- mat.mtx[
                (nrow(mat.mtx) - padSize.num + 1):nrow(mat.mtx),
            ]
        }
        mat.mtx <- rbind(mat.mtx, row.pad)
    }
    if ("left" %in% side.chr) {
        if (!is.null(value.num)) {
            col.lst <- rep(list(rep(value.num, dim(mat.mtx)[1])), padSize.num)
            col.pad <- do.call(cbind, col.lst)
        } else {
            col.pad <- mat.mtx[, padSize.num:1]
        }
        mat.mtx <- cbind(col.pad, mat.mtx)
    }
    if ("right" %in% side.chr) {
        if (!is.null(value.num)) {
            col.lst <- rep(list(rep(value.num, dim(mat.mtx)[1])), padSize.num)
            col.pad <- do.call(cbind, col.lst)
        } else {
            col.pad <- mat.mtx[,
                (ncol(mat.mtx) - padSize.num + 1):ncol(mat.mtx)
            ]
        }
        mat.mtx <- cbind(mat.mtx, col.pad)
    }
    return(mat.mtx)
}
