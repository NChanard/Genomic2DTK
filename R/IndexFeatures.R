#' Indexes GRanges on genome.
#'
#' IndexFeatures
#' @description Function that indexes a GRanges object on binned genome and constraints. Needed prior Genomic2DTK::SearchPairs() function.
#' @param gRange.gnr_lst <GRanges or GRangesList or list[GRanges]>: GRanges object, list of GRanges or GRangesList containing coordinates to index.
#' @param constraint.gnr <GRanges>: GRanges object of constraint regions. Note that bins in the same constraint region only will be paired in Genomic2DTK::SearchPairs(). If NULL chromosomes in chromSize.dtf are used as constraints (Default NULL)
#' @param chromSize.dtf <data.frame>: A data.frame containing chromosomes names and lengths in base pairs (see example).
#' @param binSize.num <integer>: Bin size in bp - corresponds to HiC matrix resolution.
#' @param variablesName.chr_vec <character> : A character vector that specify the metadata columns of GRanges on which apply the summary method if multiple ranges are indexed in the same bin.
#' @param method.chr <character>: A string defining which summary method is used on metadata columns defined in variablesName.chr_vec if multiple ranges are indexed in the same bin. Use 'mean', 'median', 'sum', 'max' or 'min'. (Default 'mean'')
#' @param cores.num <integer> : Number of cores used. (Default 1)
#' @param verbose.bln <logical>: If TRUE show the progression in console. (Default FALSE)
#' @return A GRanges object.
#' @examples
#' data(Beaf32_Peaks.gnr)
#' Beaf32_Index.gnr <- IndexFeatures(
#'     gRange.gnr_lst = list(Beaf = Beaf32_Peaks.gnr),
#'     chromSize.dtf = data.frame(
#'         seqnames = c("2L", "2R"),
#'         seqlengths = c(23513712, 25286936)
#'     ),
#'     binSize.num = 100000
#' )
#'
IndexFeatures <- function(
    gRange.gnr_lst = NULL, constraint.gnr = NULL, chromSize.dtf = NULL,
    binSize.num = NULL, method.chr = "mean", variablesName.chr_vec = NULL,
    cores.num = 1, verbose.bln = FALSE
) {
    # Constraint Informations
    if (is.null(constraint.gnr)) {
        constraint.gnr <- GenomicRanges::GRanges(
            seqnames = chromSize.dtf[, 1],
            ranges = IRanges::IRanges(
                start = rep(1, length(chromSize.dtf[, 2])),
                end = chromSize.dtf[, 2]
            ),
            strand = "*", name = chromSize.dtf[, 1]
        )
    } else {
        if (is.null(constraint.gnr$name) |
            length(which(!is.na(constraint.gnr$name))) == 0) {
            constraint.gnr$name <- paste0(
                "Constraint_",
                seq_along(constraint.gnr)
            )
        }
    }
    seqLevelsStyle.chr <- GenomeInfoDb::seqlevelsStyle(constraint.gnr)
    if (length(seqLevelsStyle.chr) > 1) {
        seqLevelsStyle.chr <- seqLevelsStyle.chr[[1]]
        GenomeInfoDb::seqlevelsStyle(constraint.gnr) <- seqLevelsStyle.chr
    }
    binnedConstraint.gnr <- BinGRanges(
        gRange.gnr = constraint.gnr,
        chromSize.dtf = chromSize.dtf,
        binSize.num = binSize.num,
        verbose.bln = verbose.bln,
        reduce.bln = FALSE,
        cores.num = cores.num
    )
    # Feature Names
    if (inherits(gRange.gnr_lst, "GRanges")) {
        gRange.gnr_lst <- list(Features = gRange.gnr_lst)
    } else if (inherits(gRange.gnr_lst, "GRangesList")) {
        gRange.gnr_lst <- as.list(gRange.gnr_lst)
    }
    if (is.null(names(gRange.gnr_lst))) {
        gRange.gnr_lst <- stats::setNames(
            gRange.gnr_lst,
            paste0("Feature_", seq_along(gRange.gnr_lst))
        )
    }
    gRangeOrder.ndx <- lapply(gRange.gnr_lst, length) |>
        unlist() |>
        order(decreasing = TRUE)
    gRange.gnr_lst <- gRange.gnr_lst[gRangeOrder.ndx]
    feature.chr_vec <- names(gRange.gnr_lst)
    # GRanges Binning
    binnedFeature.lst <- BiocParallel::bplapply(
        BPPARAM = BiocParallel::SerialParam(progressbar = verbose.bln),
        seq_along(gRange.gnr_lst),
        function(feature.ndx) {
            feature.chr <- feature.chr_vec[[feature.ndx]]
            feature.gnr <- IRanges::subsetByOverlaps(
                gRange.gnr_lst[[feature.chr]],
                constraint.gnr
            )
            GenomeInfoDb::seqlevelsStyle(feature.gnr) <- seqLevelsStyle.chr
            binnedFeature.gnr <- BinGRanges(
                gRange.gnr = feature.gnr, chromSize.dtf = chromSize.dtf,
                binSize.num = binSize.num, method.chr = method.chr,
                variablesName.chr_vec = variablesName.chr_vec,
                verbose.bln = verbose.bln, reduce.bln = TRUE,
                cores.num = cores.num
            )
            binnedFeat.tbl <- tibble::tibble(
                BinnedFeature.ndx = seq_along(binnedFeature.gnr),
                Feature.name = binnedFeature.gnr$name
            ) |>
                tidyr::unnest(cols = "Feature.name")
            binnedFeat.tbl <- dplyr::group_by(
                binnedFeat.tbl,
                Feature.name = binnedFeat.tbl$Feature.name
            )
            binnedFeat.tbl <- tidyr::nest(binnedFeat.tbl) |>
                stats::setNames(c("Feature.name", "BinnedFeature.ndx"))
            binnedConstraint.tbl <- tibble::tibble(
                BinnedConstraint.ndx = seq_along(binnedConstraint.gnr),
                Constraint.name = binnedConstraint.gnr$name
            )
            binnedConstraint.tbl <- dplyr::group_by(
                binnedConstraint.tbl,
                Constraint.name = binnedConstraint.tbl$Constraint.name
            )
            binnedConstraint.tbl <- tidyr::nest(binnedConstraint.tbl) |>
                stats::setNames(c("Constraint.name", "BinnedConstraint.ndx"))
            featConstOvlp.ovlp <- GenomicRanges::findOverlaps(
                feature.gnr,
                constraint.gnr
            )
            featConstOvlp.tbl <- tibble::tibble(
                Feature.name = feature.gnr$name[featConstOvlp.ovlp@from],
                Constraint.name = constraint.gnr$name[featConstOvlp.ovlp@to]
            ) |>
                dplyr::left_join(binnedFeat.tbl, by = "Feature.name") |>
                dplyr::select(-"Feature.name") |>
                tidyr::unnest(cols = "BinnedFeature.ndx") |>
                unique()
            featConstOvlp.tbl <- dplyr::group_by(
                featConstOvlp.tbl,
                Constraint.name = featConstOvlp.tbl$Constraint.name
            )
            featConstOvlp.tbl <- tidyr::nest(featConstOvlp.tbl) |>
                stats::setNames(c("Constraint.name", "BinnedFeature.ndx")) |>
                dplyr::left_join(binnedConstraint.tbl, by = "Constraint.name")
            multicoreParam <- MakeParallelParam(
                cores.num = cores.num,
                verbose.bln = FALSE
            )
            binnedFeature.gnr_lst <- BiocParallel::bplapply(
                BPPARAM = multicoreParam, seq_len(nrow(featConstOvlp.tbl)),
                function(row.ndx) {
                    ranges.ndx <-
                        featConstOvlp.tbl$BinnedFeature.ndx[row.ndx] |>
                        unlist(use.names = FALSE)
                    constraint.ndx <-
                        featConstOvlp.tbl$BinnedConstraint.ndx[row.ndx] |>
                        unlist(use.names = FALSE)
                    subBinnedFeature.gnr <- IRanges::subsetByOverlaps(
                        binnedFeature.gnr[ranges.ndx],
                        binnedConstraint.gnr[constraint.ndx]
                    )
                    subBinnedFeature.gnr$constraint <-
                        featConstOvlp.tbl$Constraint.name[row.ndx]
                    return(subBinnedFeature.gnr)
                }
            )
            binnedFeature.gnr <- MergeGRanges(
                binnedFeature.gnr_lst,
                sort.bln = FALSE,
                reduce.bln = FALSE
            )
            binnedFeature.gnr$bln <- TRUE
            names(S4Vectors::mcols(binnedFeature.gnr)) <- paste0(
                feature.chr,
                ".",
                names(S4Vectors::mcols(binnedFeature.gnr))
            )
            names(S4Vectors::mcols(binnedFeature.gnr))[which(
                names(S4Vectors::mcols(binnedFeature.gnr)) ==
                    paste0(feature.chr, ".bin")
            )] <- "bin"
            names(S4Vectors::mcols(binnedFeature.gnr))[which(
                names(S4Vectors::mcols(binnedFeature.gnr)) ==
                    paste0(feature.chr, ".constraint")
            )] <- "constraint"
            binnedFeature.gnr$name <- paste0(binnedFeature.gnr$bin,
                ":",
                binnedFeature.gnr$constraint
            )
            metadataBinnedFeature.dtf <- data.frame(
                S4Vectors::mcols(binnedFeature.gnr)
            )
            S4Vectors::mcols(binnedFeature.gnr) <- NULL
            return(list(
                binnedFeature.gnr = binnedFeature.gnr,
                featureMetadata.dtf = metadataBinnedFeature.dtf
            ))
        }
    )
    binnedIndex.gnr <- binnedFeature.lst |>
        lapply("[[", "binnedFeature.gnr") |>
        MergeGRanges(sort.bln = FALSE, reduce.bln = FALSE)
    featureMetadata.lst_dtf <- binnedFeature.lst |>
        lapply("[[", "featureMetadata.dtf")
    S4Vectors::mcols(binnedIndex.gnr) <- BindFillRows(featureMetadata.lst_dtf)
    ids.lst <- binnedIndex.gnr$name
    dupplicatedIds.lst <- unique(ids.lst[duplicated(ids.lst)])
    idDuplicated.ndx <- which(ids.lst %in% dupplicatedIds.lst)
    # Merge GRanges into one index and duplicated Bin handle
    if (length(idDuplicated.ndx)) {
        binnedIndexDuplicated.dtf <- data.frame(
            binnedIndex.gnr[idDuplicated.ndx]
        )
        binnedIndexDuplicated.tbl <- tibble::tibble(binnedIndexDuplicated.dtf)
        binnedIndexDuplicated.tbl <- dplyr::group_by(
            binnedIndexDuplicated.tbl,
            name = binnedIndexDuplicated.tbl$name
        )
        binnedIndexDuplicated.tbl <- tidyr::nest(binnedIndexDuplicated.tbl)
        binnedIndexNoDuplicated.dtf <- data.frame(
            binnedIndex.gnr[-idDuplicated.ndx]
        )
        binnedIndexNoDuplicated.tbl <- tibble::tibble(
            binnedIndexNoDuplicated.dtf
        )
        multicoreParam <- MakeParallelParam(
            cores.num = cores.num,
            verbose.bln = verbose.bln
        )
        binnedIndexDuplicated.lst <- BiocParallel::bplapply(
            BPPARAM = multicoreParam, seq_len(nrow(binnedIndexDuplicated.tbl)),
            function(row.ndx) {
                rowName.chr <- binnedIndexDuplicated.tbl$name[[row.ndx]]
                row <- binnedIndexDuplicated.tbl$data[[row.ndx]]
                col.lst <- lapply(
                    seq_along(row),
                    function(col.ndx) {
                        col <- dplyr::pull(row, col.ndx)
                        if (length(unique(stats::na.omit(col))) == 0) {
                            return(NA)
                        } else if (length(unique(stats::na.omit(col))) == 1) {
                            return(unique(stats::na.omit(col)))
                        } else {
                        return(list(unlist(col)))
                        }
                    }
                ) |>
                    stats::setNames(names(row))
                binnedIndexDuplicated.tbl <- tibble::as_tibble(col.lst) |>
                    tibble::add_column(name = rowName.chr)
                return(binnedIndexDuplicated.tbl)
            }
        )
        binnedIndexDuplicated.tbl <- dplyr::bind_rows(
            binnedIndexDuplicated.lst
        )
        binnedIndex.gnr <- rbind(
                binnedIndexDuplicated.tbl,
                binnedIndexNoDuplicated.tbl
            ) |>
            data.frame() |>
            methods::as("GRanges")
    }
    for (featureName.chr in feature.chr_vec) {
        colname.chr <- paste0(featureName.chr, ".bln")
        S4Vectors::mcols(binnedIndex.gnr)[
            which(is.na(S4Vectors::mcols(binnedIndex.gnr)[,colname.chr])),
            colname.chr
        ] <- FALSE
        S4Vectors::mcols(binnedIndex.gnr)[,colname.chr] <- methods::as(
            S4Vectors::mcols(binnedIndex.gnr)[,colname.chr], "Rle"
        )
    }
    columOrder.chr <- c(
            "name",
            "bin",
            "constraint",
            names(S4Vectors::mcols(binnedIndex.gnr))
        ) |>
        unique()
    S4Vectors::mcols(binnedIndex.gnr) <- as.data.frame(
            S4Vectors::mcols(binnedIndex.gnr)
        ) |>
        dplyr::select(columOrder.chr)
    return(sort(binnedIndex.gnr))
}
