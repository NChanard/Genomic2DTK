#' Compute Vanilla Count Correction.
#'
#' VCnorm
#' @description Compute Vanilla Count or Vanilla Count square root correction normalization on hic maps.
#' @param hic.cmx <contactMatrix>: The HiC maps chunk to normalize.
#' @param qtlTh.num <numerical>: the threshold quantile below which the bins will be ignored. (Default 0.15)
#' @param sqrt.bln <logical>: Whether the square root must be apply. (Default TRUE)
#' @return A matrices list.
#' @examples
#' # Data
#' data(HiC_Ctrl.cmx_lst)
#'
#' HiC_Ctrl_VC.cmx <- VCnorm(HiC_Ctrl.cmx_lst[['2L_2L']])
#' HiC_Ctrl_VC_SQRT.cmx <- VCnorm(HiC_Ctrl.cmx_lst[['2L_2L']], sqrt.bln=TRUE)

VCnorm <- function(hic.cmx = NULL, qtlTh.num = 0.15, sqrt.bln = TRUE) {
    pow.num <- ifelse(sqrt.bln, 0.5, 1)
    hic.spm <- hic.cmx@matrix
    # Removed Low counts bins
    if (qtlTh.num) {
        if (hic.cmx@metadata$symmetric) {
            rowBias.num <- Matrix::rowSums(hic.spm, na.rm = TRUE) + Matrix::colSums(hic.spm,
                na.rm = TRUE) - Matrix::diag(hic.spm)
            colBias.num <- rowBias.num
        } else {
            rowBias.num <- Matrix::rowSums(hic.spm, na.rm = TRUE)
            colBias.num <- Matrix::colSums(hic.spm, na.rm = TRUE)
        }
        row.ndx <- which(rowBias.num < stats::quantile(rowBias.num, qtlTh.num) &
            rowBias.num > 0)
        col.ndx <- which(colBias.num < stats::quantile(colBias.num, qtlTh.num) &
            colBias.num > 0)
        meltedHic.dtf <- MeltSpm(hic.spm)
        removedHic.dtf <- dplyr::filter(meltedHic.dtf, meltedHic.dtf$i %in%
            row.ndx | meltedHic.dtf$j %in% col.ndx)
        hic.cmx@metadata$removedCounts <- Matrix::sparseMatrix(i = removedHic.dtf$i,
            j = removedHic.dtf$j, x = removedHic.dtf$x, dims = dim(hic.spm))
        hic.dtf <- dplyr::filter(meltedHic.dtf, NotIn(meltedHic.dtf$i,
            row.ndx) & NotIn(meltedHic.dtf$j, col.ndx))
        hic.spm <- Matrix::sparseMatrix(i = hic.dtf$i, j = hic.dtf$j, x = hic.dtf$x,
            dims = dim(hic.spm))
    }
    hic.cmx@metadata$observed <- hic.spm@x
    # Bias computation
    coords.tbl <- MeltSpm(hic.spm)
    rowNormalizer.num <- Matrix::rowSums(hic.spm, na.rm = TRUE)
    if (hic.cmx@metadata$symmetric) {
        if (!is.na(hic.cmx@metadata$kind)) {
            rowNormalizer.num <- rowNormalizer.num + Matrix::colSums(hic.spm,
                na.rm = TRUE) - Matrix::diag(hic.spm)
        }
        colNormalizer.num <- rowNormalizer.num
    } else {
        colNormalizer.num <- Matrix::colSums(hic.spm, na.rm = TRUE)
    }
    hic.cmx@metadata$normalizer <- ((rowNormalizer.num[coords.tbl$i] *
        colNormalizer.num[coords.tbl$j])^(-pow.num)) * sum(hic.spm@x)/sum(hic.spm@x *
        ((rowNormalizer.num[coords.tbl$i] * colNormalizer.num[coords.tbl$j])^(-pow.num)))
    X.num <- hic.spm@x * (1/(rowNormalizer.num[coords.tbl$i]^pow.num)) *
        (1/(colNormalizer.num[coords.tbl$j]^pow.num))
    hic.spm@x <- X.num * sum(hic.spm@x)/sum(X.num)
    hic.cmx@matrix <- hic.spm
    hic.cmx@metadata$mtx <- "norm"
    return(hic.cmx)
}