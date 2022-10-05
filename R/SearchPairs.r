#' SearchPairs
#'
#' Compute pairs of genomic features that share a same range constraintes.
#' @param indexAnchor.gnr <GRanges>: GRanges used in the first part of the pairs.
#' @param indexBait.gnr <GRanges>: GRanges used in the second part of the pars. If NULL is replace with indexAnchor.gnr (Default NULL)
#' @param minDist.num <numeric>: the minimal distance between first and second parts of pairs. (Default NULL)
#' @param maxDist.num <numeric>: the maximal distance between first and second parts of pairs. (Default NULL)
#' @param cores.num <integer> : An integer to specify the number of cores. (Default 1)
#' @param verbose.bln <logical>: A logical value. If TRUE show the progression in console. (Default TRUE)
#' @return A GInteractions object.

SearchPairs = function(indexAnchor.gnr=NULL, indexBait.gnr=NULL, minDist.num=NULL, maxDist.num=NULL, verbose.bln=TRUE, cores.num=1){
    if(is.character(minDist.num)){
        minDist.num <- GenomicTK::GenomicSystem(minDist.num)
    }
    if(is.character(maxDist.num)){
        maxDist.num <- GenomicTK::GenomicSystem(maxDist.num)
    }
    if(is.null(indexBait.gnr)){indexBait.gnr<-indexAnchor.gnr}
    commonConstraint.lst <- intersect(indexAnchor.gnr$constraint, indexBait.gnr$constraint)
    jobLength.num <- length(commonConstraint.lst)
    if(cores.num==1){
        start.tim <- Sys.time()
        if(verbose.bln){cat("\n")}
        pairs.gni_lst <- lapply(seq_len(jobLength.num), function(constraint.ndx){
            if(verbose.bln){SuperTK::ShowLoading(start.tim, constraint.ndx,jobLength.num)}
            commonConstraint.chr <- commonConstraint.lst[[constraint.ndx]]
            subIndexAnchor.gnr <- indexAnchor.gnr[which(indexAnchor.gnr$constraint == commonConstraint.chr)]
            subIndexBait.gnr <- indexBait.gnr[which(indexBait.gnr$constraint == commonConstraint.chr)]
            pairsCombination.dtf <- expand.grid(seq_along(subIndexAnchor.gnr),seq_along(subIndexBait.gnr))
            subIndexAnchor.gnr[pairsCombination.dtf[,'Var1']]
            subIndexBait.gnr[pairsCombination.dtf[,'Var2']]
            subPairs.gni <- InteractionSet::GInteractions(subIndexAnchor.gnr[pairsCombination.dtf[,'Var1']], subIndexBait.gnr[pairsCombination.dtf[,'Var2']])
            subPairs.gni$distance <- InteractionSet::pairdist(subPairs.gni)
            if (!is.null(minDist.num)){
                subPairs.gni <- subPairs.gni[which(subPairs.gni$distance >= minDist.num)]
            }
            if (!is.null(maxDist.num)){
                subPairs.gni <- subPairs.gni[which(subPairs.gni$distance <= maxDist.num)]
            }
            return(subPairs.gni)
            }) 
        if(verbose.bln){cat("\n")}
    }else if(cores.num>=2){
        parCl <- parallel::makeCluster(cores.num, type ="FORK")
        doParallel::registerDoParallel(parCl)
        pairs.gni_lst <- parallel::parLapply(parCl,seq_len(jobLength.num), function(constraint.ndx){
            commonConstraint.chr <- commonConstraint.lst[[constraint.ndx]]
            subIndexAnchor.gnr <- indexAnchor.gnr[which(indexAnchor.gnr$constraint == commonConstraint.chr)]
            subIndexBait.gnr <- indexBait.gnr[which(indexBait.gnr$constraint == commonConstraint.chr)]
            pairsCombination.dtf <- expand.grid(seq_along(subIndexAnchor.gnr),seq_along(subIndexBait.gnr))
            subIndexAnchor.gnr[pairsCombination.dtf[,'Var1']]
            subIndexBait.gnr[pairsCombination.dtf[,'Var2']]
            subPairs.gni <- InteractionSet::GInteractions(subIndexAnchor.gnr[pairsCombination.dtf[,'Var1']], subIndexBait.gnr[pairsCombination.dtf[,'Var2']])
            subPairs.gni$distance <- InteractionSet::pairdist(subPairs.gni)
            if (!is.null(minDist.num)){
                subPairs.gni <- subPairs.gni[which(subPairs.gni$distance >= minDist.num)]
            }
            if (!is.null(maxDist.num)){
                subPairs.gni <- subPairs.gni[which(subPairs.gni$distance <= maxDist.num)]
            }
            return(subPairs.gni)
            }) 
        parallel::stopCluster(parCl)
    }
    pairs.gni <- do.call(c,pairs.gni_lst)
    pairs.gni$anchor2.constraint <- NULL
    S4Vectors::mcols(pairs.gni) %<>% as.data.frame %>% dplyr::rename(constraint=.data$anchor1.constraint)
    names(S4Vectors::mcols(pairs.gni)) <- gsub("anchor1.", "anchor.",gsub("anchor2.", "bait.",names(S4Vectors::mcols(pairs.gni))))
    S4Vectors::mcols(pairs.gni)$name <- paste0(S4Vectors::mcols(pairs.gni)$anchor.bin,"_",S4Vectors::mcols(pairs.gni)$bait.bin)
    S4Vectors::mcols(pairs.gni)$orientation <- (pairs.gni == InteractionSet::swapAnchors(pairs.gni))
    S4Vectors::mcols(pairs.gni)$submatrix.name <- S4Vectors::mcols(pairs.gni)$name
    S4Vectors::mcols(pairs.gni)$submatrix.name[!S4Vectors::mcols(pairs.gni)$orientation] <- paste0(S4Vectors::mcols(pairs.gni)$bait.bin[!S4Vectors::mcols(pairs.gni)$orientation],"_",S4Vectors::mcols(pairs.gni)$anchor.bin[!S4Vectors::mcols(pairs.gni)$orientation])
    S4Vectors::mcols(pairs.gni) %<>% as.data.frame %>% dplyr::select(.data$name,.data$constraint,.data$distance,.data$orientation,.data$submatrix.name,.data$anchor.bin,.data$anchor.name,.data$bait.bin,.data$bait.name,tidyselect::starts_with("anchor"), tidyselect::starts_with("bait")) 
    names(pairs.gni) <- S4Vectors::mcols(pairs.gni)$name
    return(pairs.gni)
}