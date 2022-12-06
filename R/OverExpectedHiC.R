#' Genomic distance bias correction.
#'
#' OverExpectedHiC
#' @description Function that normalises HiC matrices by expected values computing by genomic distance.
#' @param hic.cmx_lst <List[contactMatrix]>: The HiC maps list.
#' @param cores.num <numerical> : Number of cores to be used. (Default 1)
#' @param verbose.bln <logical>: If TRUE show the progression in console. (Default FALSE)
#' @return A matrices list.
#' @examples
#' # Note: run Genomic2DTK::BalanceHiC before OverExpectedHiC calculation.
#' data(HiC_Ctrl.cmx_lst)
#' OverExpectedHiC(HiC_Ctrl.cmx_lst)
#'
OverExpectedHiC <- function(
    hic.cmx_lst, verbose.bln = FALSE, cores.num = 1
) {
    resolution.num <- attributes(hic.cmx_lst)$resolution
    matricesNames.chr <- names(hic.cmx_lst)
    multicoreParam <- MakeParallelParam(
        cores.num = cores.num,
        verbose.bln = verbose.bln
    )
    expected.lst <- BiocParallel::bplapply(
        BPPARAM = multicoreParam, seq_along(matricesNames.chr),
        function(matrixName.ndx) {
            matrixName.chr <- matricesNames.chr[matrixName.ndx]
            if (hic.cmx_lst[[matrixName.chr]]@metadata$type == "cis") {
                hic.dtf <- MeltSpm(hic.cmx_lst[[matrixName.chr]]@matrix)
                hic.dtf <- dplyr::mutate(
                    hic.dtf,
                    distance = 1 + (hic.dtf$j - hic.dtf$i) * resolution.num) |>
                    dplyr::select(c("x", "distance"))
                expected.dtf <- dplyr::group_by(
                    hic.dtf, distance = hic.dtf$distance) |>
                    dplyr::summarise_at(
                        .vars = "x", .funs = list(expected = mean)
                    )
                expected.num <- dplyr::left_join(
                    hic.dtf,
                    expected.dtf,
                    by = "distance"
                )
            } else {
                expected.num <- list(
                    expected = sum(hic.cmx_lst[[matrixName.chr]]@matrix@x)/
                        Reduce(`*`, hic.cmx_lst[[matrixName.chr]]@matrix@Dim)
                )
            }
            expected.lst <- list(expected.num) |>
                stats::setNames(matrixName.chr)
            return(expected.lst)
        }
    ) |> do.call(what = "c")
    cisNames.chr <- NULL
    for (matrixName.chr in matricesNames.chr) {
        hic.cmx_lst[[matrixName.chr]]@metadata$expected <-
            expected.lst[[matrixName.chr]]$expected
        hic.cmx_lst[[matrixName.chr]]@matrix@x <-
            hic.cmx_lst[[matrixName.chr]]@matrix@x/
            expected.lst[[matrixName.chr]]$expected
        if (hic.cmx_lst[[matrixName.chr]]@metadata$type == "cis") {
            cisNames.chr <- c(cisNames.chr, matrixName.chr)
        }
    }
    attributes(hic.cmx_lst)$mtx <- "o/e"
    expected.dtf <- expected.lst[cisNames.chr] |>
        do.call(what = rbind)
    attributes(hic.cmx_lst)$expected <- dplyr::group_by(
        expected.dtf,
        distance = expected.dtf$distance) |>
        dplyr::summarise_at(.vars = "expected", .funs = list(expected = mean))
    return(hic.cmx_lst)
}