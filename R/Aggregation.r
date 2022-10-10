#' Aggregation
#'
#' aggregated matrices list. Could apply a differential of each paired matrices in two list before and after aggregation.
#' @param ctrlMatrices.lst <list[matrix]>: The matrices list to aggregate as control.
#' @param matrices.lst <list[matrix]>: The matrices list to aggregate.
#' @param minDist.num <numeric>: The minimal distance between anchor and bait.
#' @param maxDist.num <numeric>: The maximal distance between anchor and bait.
#' @param agg.fun <chracter or function>: The function use to aggregate each pixel in matrix. If the parameter is a character so:
#' \itemize{
#' \item "50%" or "median" apply the median
#' \item "+" or "sum" apply the sum
#' \item other (Default) apply the mean
#' }
#' @param rm0.bln <logical>: if TURE 0 are replace with NA. (Default FALSE)
#' @param diff.fun <chracter or function>: The function use to compute differential. If the parameter is character so:
#' \itemize{
#' \item "-", "substract" or "substraction" apply a substraction
#' \item "/" or "ratio" apply a ratio
#' \item "log2","log2-","log2/" or "log2ratio" apply a log2 on ratio
#' \item other (Default) apply a log2 on 1+ratio
#' }
#' @param scaleCorrection.bln <logical>: Whether a correction should be done on the median value take in ane noising area. (Default TRUE)
#' @param correctionArea.lst <list>: Nested list of indice that define a noising area fore correction. List muste contain in first an element "i" (row indices) then an element called "j" (columns indices). If NULL automatically take in upper left part of aggregated matrices. (Default NULL)
#' @param statCompare.bln <logical>: Whether a t.test must be apply to each pxl of the differential aggregated matrix.
#' @return A matrix
#' @examples
#' library(GenomicED)
#' data("submatrixRF_Ctrl.mtx_lst")
#' data("submatrixRF.mtx_lst")
#' aggreg.mtx <- Aggregation(
#'   matrices.lst = submatrixRF_Ctrl.mtx_lst, 
#'   agg.fun      = "sum",
#'   rm0.bln      = TRUE,
#'   minDist      = 9000,
#'   maxDist      = 11000
#' )
#' dim(aggreg.mtx)
#' aggreg.mtx[1:5,1:5]
#' str(attributes(aggreg.mtx),max.level = 1)
#'
#' diffAggreg.mtx <- Aggregation(
#'   ctrlMatrices.lst    = submatrixRF_Ctrl.mtx_lst,
#'   matrices.lst        = submatrixRF.mtx_lst,
#'   minDist             = 9000,
#'   maxDist             = 11000,
#'   agg.fun             = "mean",
#'   rm0.bln             = FALSE,
#'   diff.fun            = "substraction",
#'   scaleCorrection.bln = TRUE,
#'   correctionArea.lst  =  list(
#'     i = c(1:30),
#'     j = c(72:101)
#'     ),
#'   statCompare.bln = TRUE
#' )
#' dim(diffAggreg.mtx)
#' diffAggreg.mtx[1:5,1:5]
#' str(attributes(diffAggreg.mtx),max.level = 1)
Aggregation <- function(ctrlMatrices.lst=NULL, matrices.lst=NULL, minDist.num=NULL, maxDist.num=NULL, agg.fun="mean", rm0.bln=FALSE, diff.fun="substraction", scaleCorrection.bln=TRUE ,correctionArea.lst = NULL, statCompare.bln=FALSE){
    # subFunctions
        .PrepareMtxList =  function(matrices.lst, minDist.num = NULL, maxDist.num = NULL, rm0.bln=FALSE){
                interactions.gni <- attributes(matrices.lst)$interactions
            # Filter on distances
                if(!is.na(minDist.num)){
                    matrices.lst <- matrices.lst[S4Vectors::mcols(interactions.gni)$submatrix.name[which(S4Vectors::mcols(interactions.gni)$distance >= minDist.num)]]
                }
                if(!is.na(maxDist.num)){
                    matrices.lst <- matrices.lst[S4Vectors::mcols(interactions.gni)$submatrix.name[which(S4Vectors::mcols(interactions.gni)$distance <= maxDist.num)]]
                }
                matrices.lst <- matrices.lst[!is.na(names(matrices.lst))]
            # Convert sparse matrix in dense matrix and convert 0 in NA if rm0.bln is TRUE
                if(rm0.bln){
                    matrices.lst %<>% lapply(
                        function(mat.spm){
                            mat.mtx <- mat.spm %>% as.matrix
                            mat.mtx[mat.mtx==0] <- NA
                            return(mat.mtx)
                        }
                    )
                }else{
                    matrices.lst %<>% lapply(as.matrix)
                }
            #
            return(matrices.lst)
        }
    # Put list on correct variable
        if(!is.null(ctrlMatrices.lst) && is.null(matrices.lst)){
            matrices.lst <- ctrlMatrices.lst
            ctrlMatrices.lst <- NULL
        }
    # Get attributes
        matDim.num <- attributes(matrices.lst)$matriceDim
        totMtx.num  <- length(matrices.lst)
        attributes.lst <- attributes(matrices.lst)
        if("names" %in% names(attributes.lst)){attributes.lst <- attributes.lst[-which(names(attributes.lst)=="names")]}
    # Differential Function
        if(!is.function(diff.fun) && !is.null(ctrlMatrices.lst)){
            diff.fun <- dplyr::case_when(
                tolower(diff.fun) %in% c("-","substract","substraction")      ~ "function(mat.mtx,ctrl.mtx){mat.mtx - ctrl.mtx}",
                tolower(diff.fun) %in% c("/","ratio")                         ~ "function(mat.mtx,ctrl.mtx){mat.mtx / ctrl.mtx}",
                tolower(diff.fun) %in% c("log2","log2-","log2/","log2ratio")  ~ "function(mat.mtx,ctrl.mtx){log2(mat.mtx) - log2(ctrl.mtx)}",
                TRUE                                                          ~ "function(mat.mtx,ctrl.mtx){log2(mat.mtx+1) - log2(ctrl.mtx+1)}"
                )
            diff.fun <- SuperTK::WrapFunction(diff.fun)
        }
    # Aggregation Function
        if(!is.function(agg.fun)){
            agg.fun <- dplyr::case_when(
                tolower(agg.fun) %in% c("50%","median")     ~ "function(pxl){stats::median(pxl,na.rm=TRUE)}",
                tolower(agg.fun) %in% c("+","sum")          ~ "function(pxl){sum(pxl,na.rm=TRUE)}",
                TRUE                                        ~ "function(pxl){mean(pxl,na.rm=TRUE)}"
                )
            agg.fun <- SuperTK::WrapFunction(agg.fun)
        }
    # Prepare Matrix List
        if(is.null(minDist.num)){minDist.num<-NA}
        if(is.null(maxDist.num)){maxDist.num<-NA}
        matrices.lst <- .PrepareMtxList(matrices.lst=matrices.lst, minDist.num=minDist.num, maxDist.num=maxDist.num, rm0.bln=rm0.bln)
    # Aggregate
        agg.mtx <-  apply(simplify2array(matrices.lst),1:2,agg.fun)
        gc()
    # Differential Case else Return
        if(!is.null(ctrlMatrices.lst)){
            # Prepare Matrix List
                ctrlMatrices.lst <- .PrepareMtxList(matrices.lst=ctrlMatrices.lst, minDist.num=minDist.num, maxDist.num=maxDist.num, rm0.bln=rm0.bln)
            # Aggregate
                aggCtrl.mtx <- apply(simplify2array(ctrlMatrices.lst),1:2,agg.fun)
                gc()
            # Scale mat on Ctrl median
                if(scaleCorrection.bln){
                    if(is.null(correctionArea.lst) || sum(unlist(correctionArea.lst) > matDim.num)){
                        correctionArea.lst <- list(i=seq_len(round(matDim.num*0.3)), j=(matDim.num-round(matDim.num*0.3)+1):matDim.num)
                    }
                    correctionValue.num <- stats::median(aggCtrl.mtx[correctionArea.lst[[1]],correctionArea.lst[[2]]])-stats::median(agg.mtx[correctionArea.lst[[1]],correctionArea.lst[[2]]])
                    aggCorrected.mtx <- agg.mtx+correctionValue.num
                }else{
                    correctionValue.num <- NULL
                    aggCorrected.mtx <- NULL
                }
            # Stat compare
                if(statCompare.bln){
                    mtx.nlst <- matrices.lst %>% lapply(c) %>% simplify2array
                    ctrlMtx.nlst <- ctrlMatrices.lst %>% lapply(c) %>% simplify2array
                    pval.mtx <- lapply(seq_len(dim(mtx.nlst)[[1]]),function(ndx){
                        WT.vec <- mtx.nlst[ndx,]
                        KD.vec <- ctrlMtx.nlst[ndx,]
                        WT.vec <- unlist(WT.vec[!is.na(WT.vec)])
                        KD.vec <- unlist(KD.vec[!is.na(KD.vec)])
                        if(length(WT.vec)>10 & length(KD.vec)>10){
                            return(stats::t.test(WT.vec, KD.vec, var=FALSE)$p.value)
                        }else{
                            return(NA)
                        }
                    })
                    gc()
                    pval.mtx <- matrix(stats::p.adjust(c(pval.mtx), method="fdr"), nrow=matDim.num, ncol=matDim.num)
                    pval.mtx[pval.mtx>0.05] <- NA
                    pval.mtx[pval.mtx<1e-16] <- 1e-16
                    pval.mtx <- -log10(pval.mtx)
                }
            # Differential at submatrix and aggregated scale
                joblength.num <- length(ctrlMatrices.lst)
                diffmatrices.lst <- lapply(seq_len(joblength.num), function(mtx.ndx){
                    diff.mtx <- diff.fun(matrices.lst[[mtx.ndx]],ctrlMatrices.lst[[mtx.ndx]])
                    diff.mtx[is.infinite(diff.mtx)] <- NA
                    return(as.matrix(diff.mtx))
                })
                aggDelta.mtx <- diff.fun(agg.mtx,aggCtrl.mtx)
                if(!is.null(correctionArea.lst)){
                    aggCorrectedDelta.mtx <- diff.fun(aggCorrected.mtx,aggCtrl.mtx)
                }else{
                    aggCorrectedDelta.mtx <- NULL
                }
            # Aggregation of differential list
                aggDiff.mtx <- apply(simplify2array(diffmatrices.lst),1:2,agg.fun)
                gc()
            # Filter aggregated differential by pval.mtx
                if(statCompare.bln){
                    aggDiff.vec <- c(aggDiff.mtx)
                    aggDiff.vec[is.na(c(pval.mtx))] <- diff.fun(1,1)
                    aggDiffPvalFilt.mtx <- matrix(aggDiff.vec, nrow=matDim.num, ncol=matDim.num)
                }else{
                    pval.mtx <- NULL
                    aggDiffPvalFilt.mtx <- NULL
                }
            # Return
                aggDiff.mtx %>%
                    SuperTK::AddAttr(overwrite.bln=TRUE, attribute.lst=c(
                        totalMatrixNumber = totMtx.num,
                        filteredMatrixNumber = length(matrices.lst),
                        minimalDistance = minDist.num,
                        maximalDistance = maxDist.num,
                        aggregationMethod = agg.fun,
                        differentialMethod = diff.fun,
                        zeroRemoved = rm0.bln,
                        correctedFact = correctionValue.num,
                        matrices = list(list(
                            agg = agg.mtx,
                            aggCtrl = aggCtrl.mtx,
                            aggDelta = aggDelta.mtx,
                            aggCorrected = aggCorrected.mtx,
                            aggCorrectedDelta = aggCorrectedDelta.mtx,
                            pVal = pval.mtx,
                            aggDiffPvalFilt = aggDiffPvalFilt.mtx
                        )),
                        attributes.lst
                    )) %>%
                    return(.data)
        }else{
            agg.mtx %>% 
                SuperTK::AddAttr(attribute.lst=c(
                    totalMatrixNumber = totMtx.num ,
                    filteredMatrixNumber = length(matrices.lst),
                    minimalDistance = minDist.num,
                    maximalDistance = maxDist.num,
                    aggregationMethod = agg.fun,
                    zeroRemoved = rm0.bln,
                    attributes.lst
                )) %>%
                return(.data)
    }
}