
#' Merge HiC chunk.
#' 
#' JoinHiC
#' @description Create mega contactMatrix from a list of contactMatrix.
#' @param hic.cmx_lst <List[contactMatrix]>: The HiC maps list.
#' @return A contactMatrix.
#' @examples
#' data(HiC_Ctrl.cmx_lst)
#' Mega_Ctrl.cmx <- JoinHiC(HiC_Ctrl.cmx_lst)

JoinHiC <- function(hic.cmx_lst) {
    chromSize.dtf <- attributes(hic.cmx_lst)$chromSize
    mapNames <- names(hic.cmx_lst)

    mapPosition <- data.frame(apply(mapNames |>
        strsplit("_") |>
        simplify2array() |>
        t(), 1:2, function(ele) {
        which(chromSize.dtf$name == ele)
    }))
    colnames(mapPosition) <- c("i", "j")

    megaDim <- sum(chromSize.dtf$dimension)
    mega.dtf <- data.frame(i = NULL, j = NULL, x = NULL)
    for (i in seq_along(chromSize.dtf$name)) {
        addI <- 0
        if (i > 1) {
            addI <- sum(chromSize.dtf$dimension[1:i - 1])
        }
        for (j in seq(i, length(chromSize.dtf$name))) {
            addJ <- 0
            if (j > 1) {
                addJ <- sum(chromSize.dtf$dimension[1:j - 1])
            }
            if (length(which(mapPosition$i == i & mapPosition$j == j))) {
                map.dtf <- MeltSpm(hic.cmx_lst[[which(mapPosition$i ==
                  i & mapPosition$j == j)]]@matrix)
                map.dtf$i <- map.dtf$i + addI
                map.dtf$j <- map.dtf$j + addJ
                mega.dtf <- rbind(mega.dtf, map.dtf)
            }
        }
    }
    mega.spm <- Matrix::sparseMatrix(i = mega.dtf$i, j = mega.dtf$j, x = mega.dtf$x,
        dims = c(megaDim, megaDim))
    seqlengths.lst <- chromSize.dtf$length |>
        stats::setNames(chromSize.dtf$name)
    binnedGenome.gnr <- GenomicRanges::tileGenome(seqlengths.lst, tilewidth = attributes(hic.cmx_lst)$resolution,
        cut.last.tile.in.chrom = TRUE)
    megaHic.cmx <- InteractionSet::ContactMatrix(mega.spm, binnedGenome.gnr,
        binnedGenome.gnr)
    megaHic.cmx@metadata <- attributes(hic.cmx_lst)[-which(names(attributes(hic.cmx_lst)) ==
        "names")]
    megaHic.cmx@metadata$kind <- "U"
    megaHic.cmx@metadata$symmetric <- TRUE
    return(megaHic.cmx)
}