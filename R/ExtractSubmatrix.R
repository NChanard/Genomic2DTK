#' Submatrix extraction.
#'
#' ExtractSubmatrix
#' @description Extract matrices in the HiC maps list around genomic features.
#' @param feature.gn <GRanges or Pairs[GRanges] or GInteractions>: The genomic coordinates on which compute the extraction of HiC submatrix.
#' @param hic.cmx_lst <List[contactMatrix]>: The HiC maps list.
#' @param referencePoint.chr <character>: Type of extracted submatrices. "rf" for "region feature" to extract triangle-shaped matrices around regions or "pf" for "point feature" to extract square-shaped matrices around points. (Default "rf")
#' @param res.num <numeric>: The resoulution in used in hic.cmx_lst. If NULL automatically find in resolution attributes of hic.cmx_lst. (Default NULL)
#' @param matriceDim.num <numeric>: The size of matrices in output. (Default 21).
#' @param shiftFactor.num <numeric>: Only when "referencePoint.chr" is "rf". Factor defining how much of the distance between anchor and bait is extracted before and after the region (Default 1). Ex: for shiftFactor.num=2, extracted matrices will be 2*regionSize+regionSize+2*regionSize.
#' @param cores.num <integer> : An integer to specify the number of cores. (Default 1)
#' @param verbose.bln <logical>: If TRUE show the progression in console. (Default FALSE)
#' @return A matrices list.
#' @examples
#' # Data
#' data(Beaf32_Peaks.gnr)
#' data(HiC_Ctrl.cmx_lst)
#'
#' # Index Beaf32
#' Beaf32_Index.gnr <- IndexFeatures(
#'     gRange.gnr_lst = list(Beaf = Beaf32_Peaks.gnr),
#'     chromSize.dtf = data.frame(seqnames = c("2L", "2R"), seqlengths = c(23513712, 25286936)),
#'     binSize.num = 100000
#' )
#'
#' # Beaf32 <-> Beaf32 Pairing
#' Beaf_Beaf.gni <- SearchPairs(indexAnchor.gnr = Beaf32_Index.gnr)
#' Beaf_Beaf.gni <- Beaf_Beaf.gni[seq_len(2000)] # subset 2000 first for exemple
#'
#' # Matrices extractions of regions defined between Beaf32 <-> Beaf32 interactions
#' interactions_RF.mtx_lst <- ExtractSubmatrix(
#'     feature.gn = Beaf_Beaf.gni,
#'     hic.cmx_lst = HiC_Ctrl.cmx_lst,
#'     referencePoint.chr = "rf"
#' )
#'
#' # Matrices extractions center on Beaf32 <-> Beaf32 pointinteraction
#' interactions_PF.mtx_lst <- ExtractSubmatrix(
#'     feature.gn = Beaf_Beaf.gni,
#'     hic.cmx_lst = HiC_Ctrl.cmx_lst,
#'     referencePoint.chr = "pf"
#' )
#'
ExtractSubmatrix <- function(
    feature.gn = NULL, hic.cmx_lst = NULL, referencePoint.chr = "pf",
    res.num = NULL, matriceDim.num = 21, shiftFactor.num = 1, cores.num = 1,
    verbose.bln = FALSE
) {
    .GInteractionFormatting <- function(
        feature.gn, res.num
    ) {
        if (inherits(feature.gn, "GRanges")) {
            feature.gni <- InteractionSet::GInteractions(
                GenomicRanges::resize(feature.gn, 1, "start"),
                GenomicRanges::resize(feature.gn, 1, "end")
            )
        } else if (inherits(feature.gn, "Pairs") &&
            inherits(feature.gn@first, "GRanges") &&
            inherits(feature.gn@second, "GRanges")) {
            feature.gni <- InteractionSet::GInteractions(
                feature.gn@first,
                feature.gn@second
            )
        } else if (inherits(feature.gn, "GInteractions")) {
            feature.gni <- feature.gn
        }
        S4Vectors::mcols(feature.gni) <- S4Vectors::mcols(feature.gn)
        if (is.null(GenomeInfoDb::seqinfo(feature.gni))) {
            GenomeInfoDb::seqinfo(feature.gni) <-
                GenomeInfoDb::seqinfo(feature.gn)
        }
        if (is.null(S4Vectors::mcols(feature.gni)$distance)) {
            S4Vectors::mcols(feature.gni)$distance <-
            InteractionSet::pairdist(feature.gni)
        }
        if (is.null(S4Vectors::mcols(feature.gni)$orientation)) {
            S4Vectors::mcols(feature.gni)$orientation <-
                (feature.gni == InteractionSet::swapAnchors(feature.gni))
        }
        if (is.null(S4Vectors::mcols(feature.gni)$anchor.bin)) {
            S4Vectors::mcols(feature.gni)$anchor.bin <- paste0(
                GenomeInfoDb::seqnames(
                    InteractionSet::anchors(feature.gni)$first
                ), ":",
                ceiling(
                    InteractionSet::anchors(feature.gni)$first@ranges@start/
                    res.num
                )
            )
        }
        if (is.null(S4Vectors::mcols(feature.gni)$bait.bin)) {
            S4Vectors::mcols(feature.gni)$bait.bin <- paste0(
                GenomeInfoDb::seqnames(
                    InteractionSet::anchors(feature.gni)$second
                ), ":",
                ceiling(
                    InteractionSet::anchors(feature.gni)$second@ranges@start/
                    res.num
                )
            )
        }
        if (is.null(S4Vectors::mcols(feature.gni)$name)) {
            S4Vectors::mcols(feature.gni)$name <- paste0(
                S4Vectors::mcols(feature.gni)$anchor.bin, "_",
                S4Vectors::mcols(feature.gni)$bait.bin
            )
        }
        if (is.null(S4Vectors::mcols(feature.gni)$submatrix.name)) {
            S4Vectors::mcols(feature.gni)$submatrix.name <- paste0(
                S4Vectors::mcols(feature.gni)$anchor.bin, "_",
                S4Vectors::mcols(feature.gni)$bait.bin
            )
            S4Vectors::mcols(feature.gni)$submatrix.name[
                !S4Vectors::mcols(feature.gni)$orientation
                ] <- paste0(
                    S4Vectors::mcols(feature.gni)$bait.bin[
                        !S4Vectors::mcols(feature.gni)$orientation
                    ], "_",
                    S4Vectors::mcols(feature.gni)$anchor.bin[
                        !S4Vectors::mcols(feature.gni)$orientation
                    ]
                )
        }
        if (!sum(
            S4Vectors::mcols(feature.gni)$anchor.bin !=
            S4Vectors::mcols(feature.gni)$bait.bin
        )) {
            S4Vectors::mcols(feature.gni)$bait.bin <- NULL
            S4Vectors::mcols(feature.gni)$bin <-
                S4Vectors::mcols(feature.gni)$anchor.bin
            S4Vectors::mcols(feature.gni)$anchor.bin <- NULL
        }
        return(feature.gni)
    }
    # Run Check Resolution
    if (!is.null(attr(hic.cmx_lst, "resolution"))) {
        res.num <- attr(hic.cmx_lst, "resolution")
    } else if (is.character(res.num)) {
        res.num <- GenomicSystem(res.num)
    }
    # Check Dimension
    if (matriceDim.num < 5) {
        matriceDim.num <- 5
    }
    # Formatting
    feature.gn <- .GInteractionFormatting(
        feature.gn = feature.gn,
        res.num = res.num
    )
    if (!sum(feature.gn$anchor.bin != feature.gn$bait.bin)) {
        referencePoint.chr <- "pf"
    }
    # Resize Features according reference Points
    referencePoint.chr <- tolower(referencePoint.chr)
    if (referencePoint.chr == "rf") {
        cis.lgk <- ReduceRun(
            GenomeInfoDb::seqnames(InteractionSet::anchors(feature.gn)$first),
            GenomeInfoDb::seqnames(InteractionSet::anchors(feature.gn)$second),
            reduceFun.chr = "paste", sep = "_"
            ) |>
            as.character() |>
            lapply(function(combinaison) {
                split <- unlist(strsplit(combinaison, "_"))
                return(split[1] == split[2])
            }) |>
            unlist()
        feature.gn <- feature.gn[cis.lgk]
        feature.gn <- feature.gn[which(feature.gn$distance >= (3 * res.num))]
        ranges.lst_dtf <- lapply(
            c("first", "second"),
            function(anchorName.chr) {
                anchor.gnr <- InteractionSet::anchors(feature.gn)[[
                    anchorName.chr
                ]]
                anchor.dtf <- data.frame(IRanges::ranges(anchor.gnr)) |>
                    dplyr::select("start", "end")
                colnames(anchor.dtf) <- paste0(
                    anchorName.chr, ".",
                    names(anchor.dtf)
                )
                return(anchor.dtf)
            }
        )
        ranges.dtf <- do.call(cbind, ranges.lst_dtf)
        ranges.dtf <- dplyr::mutate(
            ranges.dtf,
            start = pmin(ranges.dtf$first.start, ranges.dtf$second.start)
        )
        ranges.dtf <- dplyr::mutate(
            ranges.dtf,
            end = pmax(ranges.dtf$first.end, ranges.dtf$second.end)
        )
        feature.gnr <- GenomicRanges::GRanges(
            seqnames = GenomeInfoDb::seqnames(
                InteractionSet::anchors(feature.gn)$first
            ),
            ranges = IRanges::IRanges(
                start = ranges.dtf$start,
                end = ranges.dtf$end
            ),
            seqlengths = GenomeInfoDb::seqlengths(feature.gn),
            seqinfo = GenomeInfoDb::seqinfo(feature.gn)
        )
        featureResize.gnr <- GenomicRanges::resize(
            feature.gnr,
            width = feature.gnr@ranges@width +
                feature.gn$distance * shiftFactor.num * 2,
            fix = "center"
        )
        featureResize.gni <- InteractionSet::GInteractions(
            featureResize.gnr,
            featureResize.gnr
        )
        S4Vectors::mcols(featureResize.gni) <- S4Vectors::mcols(feature.gn)
    } else if (referencePoint.chr == "pf") {
        featureResize.gni <- suppressWarnings(
            GenomicRanges::resize(
                feature.gn, width = res.num * (matriceDim.num - 1) + 1,
                fix = "center"
            )
        )
    }
    # Filt Out Of Bound
    featureFilt.gni <- featureResize.gni[which(
        1L <=
        data.frame(
            InteractionSet::anchors(featureResize.gni)$first@ranges
        )[,"start"] &
        data.frame(
            InteractionSet::anchors(featureResize.gni)$first@ranges
        )[,"end"] <=
        SeqEnds(InteractionSet::anchors(featureResize.gni)$first) &
        1L <=
        data.frame(
            InteractionSet::anchors(featureResize.gni)$second
        )[,"start"] &
        data.frame(
            InteractionSet::anchors(featureResize.gni)$second
        )[,"end"] <=
        SeqEnds(InteractionSet::anchors(featureResize.gni)$second)
    )]
    # Filt Duplicated Submatrix before extraction
    featureNoDup.gni <-
        featureFilt.gni[!duplicated(featureFilt.gni$submatrix.name)]
    # Order according Chromosomes combinaison
    chromosomesCombinaison.rle <- ReduceRun(
        GenomeInfoDb::seqnames(
            InteractionSet::anchors(featureNoDup.gni)$first
        ),
        GenomeInfoDb::seqnames(
            InteractionSet::anchors(featureNoDup.gni)$second
        ),
        reduceFun.chr = "paste", sep = "_"
    )
    order.num <- unlist(
        lapply(
            unique(S4Vectors::runValue(chromosomesCombinaison.rle)),
            function(comb) {
                which(as.character(chromosomesCombinaison.rle) == comb)
            }
        )
    )
    featureNoDup.gni <- featureNoDup.gni[order.num]
    chromosomesCombinaison.rle <- chromosomesCombinaison.rle[order.num]
    matAnchors.gnr_lst <- lapply(hic.cmx_lst, InteractionSet::anchors)
    # Keep combinaison that are present in hic.cmx_lst
    present.ndx <- which(
        as.vector(chromosomesCombinaison.rle) %in% names(matAnchors.gnr_lst) |
            as.vector(chromosomesCombinaison.rle) %in% unlist(
                lapply(
                    names(matAnchors.gnr_lst),
                    function(name.chr) {
                        strsplit(name.chr, "_") |>
                        unlist() |>
                        rev() |>
                        paste(collapse = "_")
                    }
                )
            )
    )
    chromosomesCombinaison.rle <- chromosomesCombinaison.rle[present.ndx]
    featureNoDup.gni <- featureNoDup.gni[present.ndx]
    # Separate anchors
    anchors.gnr <- InteractionSet::anchors(featureNoDup.gni)$first
    baits.gnr <- InteractionSet::anchors(featureNoDup.gni)$second
    # Extraction
    submatrix.spm_lst <- BiocParallel::bplapply(
        BPPARAM = BiocParallel::SerialParam(progressbar = verbose.bln),
        seq_along(S4Vectors::runValue(chromosomesCombinaison.rle)),
        function(combinaison.ndx) {
            combinaisonName.chr <-
                S4Vectors::runValue(
                    chromosomesCombinaison.rle
                )[[combinaison.ndx]]
            combinaisonStart.ndx <-
                cumsum(
                    c(1, S4Vectors::runLength(chromosomesCombinaison.rle))
                )[[combinaison.ndx]]
            combinaisonEnd.ndx <-
                cumsum(
                    S4Vectors::runLength(chromosomesCombinaison.rle)
                )[[combinaison.ndx]]
            if (combinaisonName.chr %in% names(matAnchors.gnr_lst)) {
                mat.ndx <- which(names(hic.cmx_lst) == combinaisonName.chr)
                ovl_row <- data.frame(GenomicRanges::findOverlaps(
                    anchors.gnr[combinaisonStart.ndx:combinaisonEnd.ndx],
                    matAnchors.gnr_lst[[mat.ndx]]$row
                ))
                ovl_row <- dplyr::group_by(
                    ovl_row,
                    queryHits = ovl_row$queryHits
                )
                ovl_row <- tidyr::nest(ovl_row)
                ovl_col <- data.frame(GenomicRanges::findOverlaps(
                    baits.gnr[combinaisonStart.ndx:combinaisonEnd.ndx],
                    matAnchors.gnr_lst[[mat.ndx]]$col
                ))
                ovl_col <- dplyr::group_by(
                    ovl_col,
                    queryHits = ovl_col$queryHits
                )
                ovl_col <- tidyr::nest(ovl_col)
            } else if ({
                combinaisonName.chr |>
                strsplit("_") |>
                unlist() |>
                rev() |>
                paste(collapse = "_")} %in% names(matAnchors.gnr_lst)) {
                    mat.ndx <- which(
                        names(hic.cmx_lst) == {
                            combinaisonName.chr |>
                            strsplit("_") |>
                            unlist() |>
                            rev() |>
                            paste(collapse = "_")
                        }
                    )
                    ovl_row <- data.frame(GenomicRanges::findOverlaps(
                        baits.gnr[combinaisonStart.ndx:combinaisonEnd.ndx],
                        matAnchors.gnr_lst[[mat.ndx]]$row
                    ))
                    ovl_row <- dplyr::group_by(
                        ovl_row,
                        queryHits = ovl_row$queryHits
                    )
                    ovl_row <- tidyr::nest(ovl_row)
                    ovl_col <- data.frame(GenomicRanges::findOverlaps(
                        anchors.gnr[combinaisonStart.ndx:combinaisonEnd.ndx],
                        matAnchors.gnr_lst[[mat.ndx]]$col
                    ))
                    ovl_col <- dplyr::group_by(
                        ovl_col,
                        queryHits = ovl_col$queryHits
                    )
                    ovl_col <- tidyr::nest(ovl_col)
            }
            multicoreParam <- MakeParallelParam(
                cores.num = cores.num,
                verbose.bln = FALSE
            )
            tempSubmatrix.spm_lst <- BiocParallel::bplapply(
                BPPARAM = multicoreParam,
                seq_along(combinaisonStart.ndx:combinaisonEnd.ndx),
                function(range.ndx) {
                    row.ndx <- unlist(
                        ovl_row[[range.ndx, "data"]],
                        use.names = FALSE
                    )
                    col.ndx <- unlist(
                        ovl_col[[range.ndx, "data"]],
                        use.names = FALSE
                    )
                    if (S4Vectors::metadata(hic.cmx_lst[[mat.ndx]])$type ==
                        "cis"
                    ) {
                        gap.num <-
                            stats::median(col.ndx) - stats::median(row.ndx)
                    } else {
                        gap.num <- Inf
                    }
                    if (gap.num < 0) {
                        mat.spm <-
                            hic.cmx_lst[[mat.ndx]][col.ndx, row.ndx]@matrix
                    } else {
                        mat.spm <-
                            hic.cmx_lst[[mat.ndx]][row.ndx, col.ndx]@matrix
                    }
                    if (dim(mat.spm)[1] != matriceDim.num) {
                        mat.spm <- ResizeMatrix(
                            matrice.mtx = mat.spm,
                            newDim.num = c(matriceDim.num, matriceDim.num)
                        )
                    }
                    if (abs(gap.num) < matriceDim.num) {
                        if (abs(gap.num) > 0) {
                            mat.spm <- PadMtx(
                                mat.mtx = mat.spm, padSize.num = abs(gap.num),
                                value.num = 0, side.chr = c("left", "bot")
                            )
                        }
                        mat.spm[lower.tri(mat.spm)] <- NA
                        if (abs(gap.num) > 0) {
                            mat.spm <- mat.spm[
                                seq_len(matriceDim.num),
                                (abs(gap.num)+1):(matriceDim.num+abs(gap.num))
                            ]
                        }
                    }
                    return(as.matrix(mat.spm))
                }
            )
            return(tempSubmatrix.spm_lst)
        }
    ) |>
        do.call(what = c) |>
        stats::setNames(featureNoDup.gni$submatrix.name)
    submatrix.spm_lst <- submatrix.spm_lst[featureFilt.gni$submatrix.name]
    interactions.ndx <- seq_along(feature.gn$name) |>
        stats::setNames(feature.gn$name)
    interactions.ndx <- interactions.ndx[featureFilt.gni$name]
    attributes(submatrix.spm_lst)$interactions <- feature.gn[interactions.ndx]
    attributes(submatrix.spm_lst)$resolution <- res.num
    attributes(submatrix.spm_lst)$referencePoint <- referencePoint.chr
    attributes(submatrix.spm_lst)$matriceDim <- matriceDim.num
    if (referencePoint.chr == "rf") {
        attributes(submatrix.spm_lst)$shiftFactor <- shiftFactor.num
    }
    return(submatrix.spm_lst)
}
