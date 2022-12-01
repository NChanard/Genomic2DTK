#' Change values of HiC map.
#'
#' SwitchMatrix
#' @description Change values in matrix with observed, balanced, observed/expected or expected values according to what are be done in hic.
#' @param hic.cmx_lst <List[contactMatrix]>: The HiC maps list.
#' @param matrixKind.chr <character>: The kind of matrix you want.
#' @return A contactMatrix list.
#' @examples
#' # Data
#' data(HiC_Ctrl.cmx_lst)
#'
#' # Preprocess HiC
#' HiC.cmx_lst <- HiC_Ctrl.cmx_lst |>
#'     NormalizeHiC(
#'         interaction.type="cis",
#'         method.chr="ICE") |> 
#'     ExpectedHiC()
#' 
#' # Switch values in matrix
#' HiC_Ctrl_Obs.cmx_lst <- SwitchMatrix(HiC.cmx_lst, matrixKind.chr="obs")
#' HiC_Ctrl_Norm.cmx_lst <- SwitchMatrix(HiC.cmx_lst, matrixKind.chr="norm")
#' HiC_Ctrl_oe.cmx_lst <- SwitchMatrix(HiC.cmx_lst, matrixKind.chr="o/e")
#' HiC_Ctrl_exp.cmx_lst <- SwitchMatrix(HiC.cmx_lst, matrixKind.chr="exp")


SwitchMatrix <- function(hic.cmx_lst, matrixKind.chr = c("obs", "norm",
    "o/e", "exp")) {
    if (!(matrixKind.chr %in% c("obs", "norm", "o/e", "exp"))) {
        stop("ERROR: matrixKind.chr must be one of \"obs\", \"norm\", \"o/e\", \"exp\".")
    }
    if (attributes(hic.cmx_lst)$mtx != matrixKind.chr) {
        lapply(names(hic.cmx_lst), function(hicName.chr) {
            hic.cmx_lst[[hicName.chr]]@matrix@x <<- dplyr::case_when(matrixKind.chr ==
                "obs" ~ (hic.cmx_lst[[hicName.chr]]@metadata$observed),
                matrixKind.chr == "norm" ~ (hic.cmx_lst[[hicName.chr]]@metadata$observed *
                  hic.cmx_lst[[hicName.chr]]@metadata$normalizer), matrixKind.chr ==
                  "o/e" ~ (hic.cmx_lst[[hicName.chr]]@metadata$observed *
                  hic.cmx_lst[[hicName.chr]]@metadata$normalizer/hic.cmx_lst[[hicName.chr]]@metadata$expected),
                matrixKind.chr == "exp" ~ (hic.cmx_lst[[hicName.chr]]@metadata$expected))
        })
        attributes(hic.cmx_lst)$mtx <- matrixKind.chr
        return(hic.cmx_lst)
    } else {
        message(paste0("\nhic.cmx_lst is already ", matrixKind.chr, ".\n"))
    }
}