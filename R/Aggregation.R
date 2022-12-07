#' Aggregation of matrices list.
#'
#' Aggregation
#' @description Aggregates all the matrices of a list (or two lists in case of differential aggregation) into a single matrix. This function allows to apply different aggregation (average, sum, ...), transformation (rank, percentage, ...) and differential (subtraction, ratio, ...) functions.
#' @param ctrlMatrices.lst <list[matrix]>: The matrices list to aggregate as control.
#' @param matrices.lst <list[matrix]>: The matrices list to aggregate.
#' @param minDist.num <numeric>: The minimal distance between anchor and bait.
#' @param maxDist.num <numeric>: The maximal distance between anchor and bait.
#' @param agg.fun <function or chracter>: The function use to aggregate each pixel in matrix. If the parameter is a character so:
#' \itemize{
#' \item "50%" or "median" apply the median
#' \item "+" or "sum" apply the sum
#' \item other (Default) apply the mean
#' }
#' @param rm0.bln <logical>: Whether 0 are replace with NA. (Default FALSE)
#' @param diff.fun <function or chracter>: The function use to compute differential. If the parameter is character so:
#' \itemize{
#' \item "-", "substract" or "substraction" apply a substraction
#' \item "/" or "ratio" apply a ratio
#' \item "log2","log2-","log2/" or "log2ratio" apply a log2 on ratio
#' \item other (Default) apply a log2 on 1+ratio
#' }
#' @param trans.fun <function or chracter>: The function use to transforme or scale values in each submatrix before aggregation. If the parameter is character so:
#' \itemize{
#' \item "quantile" or "qtl" apply function dplyr::ntile(x,500)
#' \item "percentile" or "prct" apply percentile.
#' \item "rank" apply a ranking.
#' \item "zscore" apply a scaling.
#' \item "minmax" apply a Genomic2DTK::MinMaxScale.
#' \item "mu" apply a Genomic2DTK::MeanScale.
#' \item other or NULL don't apply transformation (Default).
#' }
#' @param scaleCorrection.bln <logical>: Whether a correction should be done on the median value take in ane noising area. (Default TRUE)
#' @param correctionArea.lst <list>: Nested list of indice that define a noising area fore correction. List muste contain in first an element "i" (row indices) then an element called "j" (columns indices). If NULL automatically take in upper left part of aggregated matrices. (Default NULL)
#' @param statCompare.bln <logical>: Whether a t.test must be apply to each pxl of the differential aggregated matrix.
#' @return A matrix
#' @examples
#' # Data
#' data(Beaf32_Peaks.gnr)
#' data(HiC_Ctrl.cmx_lst)
#' data(HiC_HS.cmx_lst)
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
#' # Matrices extractions center on Beaf32 <-> Beaf32 point interaction
#' interactions_Ctrl.mtx_lst <- ExtractSubmatrix(
#'     feature.gn = Beaf_Beaf.gni,
#'     hic.cmx_lst = HiC_Ctrl.cmx_lst,
#'     referencePoint.chr = "pf"
#' )
#' interactions_HS.mtx_lst <- ExtractSubmatrix(
#'     feature.gn = Beaf_Beaf.gni,
#'     hic.cmx_lst = HiC_HS.cmx_lst,
#'     referencePoint.chr = "pf"
#' )
#'
#' # Aggregate matrices in one matrix
#' aggreg.mtx <- Aggregation(interactions_Ctrl.mtx_lst)
#'
#' # Differential Aggregation
#' aggregDiff.mtx <- Aggregation(
#'     ctrlMatrices.lst = interactions_Ctrl.mtx_lst,
#'     matrices.lst = interactions_HS.mtx_lst
#' )
#'
Aggregation <- function(
    ctrlMatrices.lst = NULL, matrices.lst = NULL, minDist.num = NULL,
    maxDist.num = NULL, trans.fun = NULL, agg.fun = "mean", rm0.bln = FALSE,
    diff.fun = "substraction", scaleCorrection.bln = FALSE,
    correctionArea.lst = NULL, statCompare.bln = FALSE
) {
    # subFunctions
    .PrepareMtxList <- function(
        matrices.lst, minDist.num = NULL, maxDist.num = NULL, rm0.bln = FALSE,
        trans.fun = NULL
    ) {
        interactions.gni <- attributes(matrices.lst)$interactions
        # Filter on distances
        if (!is.na(minDist.num)) {
            if (is.character(minDist.num)) {
                minDist.num <- GenomicSystem(minDist.num)
            }
            matrices.lst <- matrices.lst[
                S4Vectors::mcols(interactions.gni)$submatrix.name[which(
                    S4Vectors::mcols(interactions.gni)$distance >= minDist.num
                )]
            ]
        }
        if (!is.na(maxDist.num)) {
            if (is.character(maxDist.num)) {
                maxDist.num <- GenomicSystem(maxDist.num)
            }
            matrices.lst <- matrices.lst[
                S4Vectors::mcols(interactions.gni)$submatrix.name[which(
                    S4Vectors::mcols(interactions.gni)$distance <= maxDist.num
                )]
            ]
        }
        matrices.lst <- matrices.lst[!is.na(names(matrices.lst))]
        # Convert sparse matrix in dense matrix and convert 0 in NA
        # if rm0.bln is TRUE
        matrices.lst <- lapply(
            matrices.lst, function(mat.spm) {
                mat.mtx <- as.matrix(mat.spm)
                if (rm0.bln) {
                    mat.mtx[mat.mtx == 0] <- NA
                }
                if (!is.null(trans.fun)) {
                    mat.mtx <- trans.fun(mat.mtx)
                }
                return(mat.mtx)
            }
        )
        return(matrices.lst)
    }
    # Put list on correct variable
    if (!is.null(ctrlMatrices.lst) && is.null(matrices.lst)) {
        matrices.lst <- ctrlMatrices.lst
        ctrlMatrices.lst <- NULL
    }
    # Get attributes
    matDim.num <- attributes(matrices.lst)$matriceDim
    totMtx.num <- length(matrices.lst)
    attributes.lst <- attributes(matrices.lst)
    if ("names" %in% names(attributes.lst)) {
        attributes.lst <- attributes.lst[-which(
            names(attributes.lst) == "names"
        )]
    }
    # Differential Function
    if (!is.function(diff.fun) &&
        !is.null(ctrlMatrices.lst)) {
        diff.fun <- dplyr::case_when(
            tolower(diff.fun) %in% c("-", "substract", "substraction") ~
                "function(mat.mtx,ctrl.mtx){mat.mtx - ctrl.mtx}",
            tolower(diff.fun) %in% c("/", "ratio") ~
                "function(mat.mtx,ctrl.mtx){mat.mtx / ctrl.mtx}",
            tolower(diff.fun) %in% c("log2", "log2-", "log2/", "log2ratio") ~
                "function(mat.mtx,ctrl.mtx){log2(mat.mtx)-log2(ctrl.mtx)}",
            TRUE ~
                "function(mat.mtx,ctrl.mtx){log2(mat.mtx+1)-log2(ctrl.mtx+1)}"
        )
        diff.fun <- WrapFunction(diff.fun)
    }
    # Aggregation Function
    if (!is.function(agg.fun)) {
        agg.fun <- dplyr::case_when(
            tolower(agg.fun) %in% c("50%", "median") ~
                "function(pxl){stats::median(pxl,na.rm=TRUE)}",
            tolower(agg.fun) %in% c("+", "sum") ~
                "function(pxl){sum(pxl,na.rm=TRUE)}",
            TRUE ~
                "function(pxl){mean(pxl,na.rm=TRUE,trim=0.01)}"
        )
        agg.fun <- WrapFunction(agg.fun)
    }
    # Transformation Function
    if (!is.function(trans.fun) &
        !is.null(trans.fun)) {
        trans.fun <- dplyr::case_when(
            tolower(trans.fun) %in% c("quantile", "qtl") ~
                "function(mat.mtx){
                    matrix(
                        dplyr::ntile(mat.mtx,500),
                        dim(mat.mtx)[[1]],
                        dim(mat.mtx)[[2]]
                    )
                }",
            tolower(trans.fun) %in% c("percentile", "prct") ~
                "function(mat.mtx){
                    matrix(
                        rank(mat.mtx,na.last='keep')/
                        length(mat.mtx[!is.na(mat.mtx)]),
                        dim(mat.mtx)[[1]],
                        dim(mat.mtx)[[2]]
                    )
                }",
            tolower(trans.fun) %in% c("rank") ~
                "function(mat.mtx){
                    matrix(
                        rank(mat.mtx,na.last='keep'),
                        dim(mat.mtx)[[1]],
                        dim(mat.mtx)[[2]]
                    )
                }",
            tolower(trans.fun) %in% c("zscore") ~
                "function(mat.mtx){
                    matrix(
                        scale(c(mat.mtx)),
                        dim(mat.mtx)[[1]],
                        dim(mat.mtx)[[2]]
                    )
                }",
            tolower(trans.fun) %in% c("minmax") ~
                "function(mat.mtx){
                    matrix(
                        MinMaxScale(c(mat.mtx)),
                        dim(mat.mtx)[[1]],
                        dim(mat.mtx)[[2]]
                    )
                }",
            tolower(trans.fun) %in% c("mu") ~
                "function(mat.mtx){
                    matrix(
                        MeanScale(c(mat.mtx)),
                        dim(mat.mtx)[[1]],
                        dim(mat.mtx)[[2]]
                    )
                }",
            TRUE ~
                "NULL"
        )
        trans.fun <- WrapFunction(trans.fun)
    }
    # Prepare Matrix List
    if (is.null(minDist.num)) {
        minDist.num <- NA
    }
    if (is.null(maxDist.num)) {
        maxDist.num <- NA
    }
    matrices.lst <- .PrepareMtxList(
        matrices.lst = matrices.lst, minDist.num = minDist.num,
        maxDist.num = maxDist.num, rm0.bln = rm0.bln, trans.fun = trans.fun
    )
    # Aggregate
    agg.mtx <- apply(
        simplify2array(matrices.lst),
        seq_len(2),
        agg.fun
    )
    gc()
    # Differential Case else Return
    if (!is.null(ctrlMatrices.lst)) {
        # Prepare Matrix List
        ctrlMatrices.lst <- .PrepareMtxList(
            matrices.lst = ctrlMatrices.lst, minDist.num = minDist.num,
            maxDist.num = maxDist.num, rm0.bln = rm0.bln, trans.fun = trans.fun
        )
        # Aggregate
        aggCtrl.mtx <- apply(
            simplify2array(ctrlMatrices.lst),
            seq_len(2),
            agg.fun
        )
        gc()
        # Scale mat on Ctrl median
        if (scaleCorrection.bln) {
            if (is.null(correctionArea.lst) ||
            sum(unlist(correctionArea.lst) > matDim.num)) {
                correctionArea.lst <- list(
                    i = seq_len(round(matDim.num * 0.3)),
                    j = (matDim.num - round(matDim.num * 0.3) + 1):matDim.num
                )
            }
            correctionValue.num <- stats::median(
                    aggCtrl.mtx[correctionArea.lst[[1]],
                    correctionArea.lst[[2]]]) -
                stats::median(
                    agg.mtx[correctionArea.lst[[1]],
                    correctionArea.lst[[2]]])
            aggCorrected.mtx <- agg.mtx + correctionValue.num
        } else {
            correctionValue.num <- NULL
            aggCorrected.mtx <- NULL
        }
        # Stat compare
        if (statCompare.bln) {
            mtx.nlst <- simplify2array(lapply(matrices.lst, c))
            ctrlMtx.nlst <- simplify2array(lapply(ctrlMatrices.lst, c))
            pval.mtx <- lapply(
                seq_len(dim(mtx.nlst)[[1]]),
                function(ndx) {
                    WT.vec <- mtx.nlst[ndx, ]
                    KD.vec <- ctrlMtx.nlst[ndx, ]
                    WT.vec <- unlist(WT.vec[!is.na(WT.vec)])
                    KD.vec <- unlist(KD.vec[!is.na(KD.vec)])
                    if (length(WT.vec) > 10 & length(KD.vec) > 10) {
                        return(stats::t.test(WT.vec, KD.vec, var = FALSE)$p.value)
                    } else {
                        return(NA)
                    }
                }
            )
            gc()
            pval.mtx <- matrix(
                stats::p.adjust(
                    c(pval.mtx),
                    method = "fdr"
                ),
                nrow = matDim.num,
                ncol = matDim.num
            )
            pval.mtx[pval.mtx > 0.05] <- NA
            pval.mtx[pval.mtx < 1e-16] <- 1e-16
            pval.mtx <- -log10(pval.mtx)
        }
        # Differential at submatrix and aggregated scale
        diffmatrices.lst <- lapply(
            seq_along(ctrlMatrices.lst),
            function(mtx.ndx) {
                diff.mtx <- diff.fun(
                    matrices.lst[[mtx.ndx]],
                    ctrlMatrices.lst[[mtx.ndx]])
                diff.mtx[is.infinite(diff.mtx)] <- NA
                return(as.matrix(diff.mtx))
            }
        )
        aggDelta.mtx <- diff.fun(agg.mtx, aggCtrl.mtx)
        if (!is.null(correctionArea.lst)) {
            aggCorrectedDelta.mtx <- diff.fun(aggCorrected.mtx, aggCtrl.mtx)
        } else {
            aggCorrectedDelta.mtx <- NULL
        }
        # Aggregation of differential list
        aggDiff.mtx <- apply(
            simplify2array(diffmatrices.lst),
            seq_len(2),
            agg.fun
        )
        gc()
        # Filter aggregated differential by pval.mtx
        if (statCompare.bln) {
            aggDiff.vec <- c(aggDiff.mtx)
            aggDiff.vec[is.na(c(pval.mtx))] <- diff.fun(1, 1)
            aggDiffPvalFilt.mtx <- matrix(
                aggDiff.vec,
                nrow = matDim.num,
                ncol = matDim.num)
        } else {
            pval.mtx <- NULL
            aggDiffPvalFilt.mtx <- NULL
        }
        # Return
        aggDiff.mtx <- AddAttr(
            var.any = aggDiff.mtx,
            overwrite.bln = TRUE,
            attribute.lst = c(
                totalMatrixNumber = totMtx.num,
                filteredMatrixNumber = length(matrices.lst),
                minimalDistance = minDist.num, maximalDistance = maxDist.num,
                transformationMethod = trans.fun, aggregationMethod = agg.fun,
                differentialMethod = diff.fun, zeroRemoved = rm0.bln,
                correctedFact = correctionValue.num,
                correctionArea = correctionArea.lst,
                matrices = list(list(
                    agg = agg.mtx, aggCtrl = aggCtrl.mtx, aggDelta = aggDelta.mtx,
                    aggCorrected = aggCorrected.mtx,
                    aggCorrectedDelta = aggCorrectedDelta.mtx, pVal = pval.mtx,
                    aggDiffPvalFilt = aggDiffPvalFilt.mtx
                )),
                attributes.lst
            )
        )
        return(aggDiff.mtx)
    } else {
        agg.mtx <- AddAttr(
            var.any = agg.mtx,
            attribute.lst = c(
                totalMatrixNumber = totMtx.num,
                filteredMatrixNumber = length(matrices.lst),
                minimalDistance = minDist.num, maximalDistance = maxDist.num,
                transformationMethod = trans.fun, aggregationMethod = agg.fun,
                zeroRemoved = rm0.bln, attributes.lst
            )
        )
        return(agg.mtx)
    }
}
