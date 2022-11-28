#' Creates pairs of coordinates from indexed anchor and bait genomic coordinates according to distance constraints.
#' 
#' SearchPairs
#' @description Creates pairs of coordinates from indexed anchor and bait genomic coordinates according to distance constraints.
#' @param indexAnchor.gnr <GRanges>: A first indexed GRanges object used as pairs anchor (must be indexed using IndexFeatures()).
#' @param indexBait.gnr <GRanges>: A second indexed GRanges object used as pairs bait (must be indexed using IndexFeatures()). If NULL, indexAnchor.gnr is used instead (Default NULL)
#' @param minDist.num <numeric>: Minimal distance between anchors and baits. (Default NULL)
#' @param maxDist.num <numeric>: Maximal distance between anchors and baits. (Default NULL)
#' @param cores.num <integer> : Number of cores to use. (Default 1)
#' @param verbose.bln <logical>: If TRUE show the progression in console. (Default TRUE)
#' @return A GInteractions object.
#' @examples
#' library(GenomicED)
#' 
#' data(anchors_Index.gnr)
#' anchors_Index.gnr[1]
#' 
#' data(baits_Index.gnr)
#' baits_Index.gnr[1]
#' 
#' \dontrun{
#' interactions.gni <- SearchPairs(
#'     indexAnchor.gnr = anchors_Index.gnr,
#'     indexBait.gnr   = baits_Index.gnr,
#'     minDist.num     = 9000, 
#'     maxDist.num     = 11000,
#'     cores.num       = 1,
#'     verbose.bln     = FALSE
#'     )
#' 
#' }

SearchPairs = function(indexAnchor.gnr=NULL, indexBait.gnr=NULL, minDist.num=NULL, maxDist.num=NULL, verbose.bln=TRUE, cores.num=1){
    if(is.character(minDist.num)){
        minDist.num <- GenomicSystem(minDist.num)
    }
    if(is.character(maxDist.num)){
        maxDist.num <- GenomicSystem(maxDist.num)
    }
    if(is.null(indexBait.gnr)){indexBait.gnr<-indexAnchor.gnr}
    commonConstraint.lst <- intersect(indexAnchor.gnr$constraint, indexBait.gnr$constraint)
    jobLength.num <- length(commonConstraint.lst)
    # if(cores.num==1){
    #     start.tim <- Sys.time()
    #     if(verbose.bln){cat("\n")}
    #     pairs.gni_lst <- lapply(seq_len(jobLength.num), function(constraint.ndx){
    #         if(verbose.bln){ShowLoading(start.tim, constraint.ndx,jobLength.num)}
    #         commonConstraint.chr <- commonConstraint.lst[[constraint.ndx]]
    #         subIndexAnchor.gnr <- indexAnchor.gnr[which(indexAnchor.gnr$constraint == commonConstraint.chr)]
    #         subIndexBait.gnr <- indexBait.gnr[which(indexBait.gnr$constraint == commonConstraint.chr)]
    #         pairsCombination.dtf <- expand.grid(seq_along(subIndexAnchor.gnr),seq_along(subIndexBait.gnr))
    #         subIndexAnchor.gnr[pairsCombination.dtf[,'Var1']]
    #         subIndexBait.gnr[pairsCombination.dtf[,'Var2']]
    #         subPairs.gni <- InteractionSet::GInteractions(subIndexAnchor.gnr[pairsCombination.dtf[,'Var1']], subIndexBait.gnr[pairsCombination.dtf[,'Var2']])
    #         subPairs.gni$distance <- InteractionSet::pairdist(subPairs.gni)
    #         if (!is.null(minDist.num)){
    #             subPairs.gni <- subPairs.gni[which(subPairs.gni$distance >= minDist.num)]
    #         }
    #         if (!is.null(maxDist.num)){
    #             subPairs.gni <- subPairs.gni[which(subPairs.gni$distance <= maxDist.num)]
    #         }
    #         return(subPairs.gni)
    #         }) 
    #     if(verbose.bln){cat("\n")}
    # }else if(cores.num>=2){
        multicoreParam <- makeParallelParam(cores.num = cores.num, verbose.bln = verbose.bln)
        pairs.gni_lst <- BiocParallel::bplapply(BPPARAM = multicoreParam,seq_len(jobLength.num), function(constraint.ndx){
            commonConstraint.chr <- commonConstraint.lst[[constraint.ndx]]
            subIndexAnchor.gnr <- indexAnchor.gnr[which(indexAnchor.gnr@elementMetadata$constraint == commonConstraint.chr)]
            subIndexBait.gnr <- indexBait.gnr[which(indexBait.gnr@elementMetadata$constraint == commonConstraint.chr)]
            pairsCombination.dtf <- expand.grid(seq_along(subIndexAnchor.gnr),seq_along(subIndexBait.gnr))
            subIndexAnchor.gnr[pairsCombination.dtf[,'Var1']]
            subIndexBait.gnr[pairsCombination.dtf[,'Var2']]
            subPairs.gni <- InteractionSet::GInteractions(subIndexAnchor.gnr[pairsCombination.dtf[,'Var1']], subIndexBait.gnr[pairsCombination.dtf[,'Var2']])
            subPairs.gni$distance <- InteractionSet::pairdist(subPairs.gni)
            if (!is.null(minDist.num)){
                subPairs.gni <- subPairs.gni[which(subPairs.gni@elementMetadata$distance >= minDist.num)]
            }
            if (!is.null(maxDist.num)){
                subPairs.gni <- subPairs.gni[which(subPairs.gni@elementMetadata$distance <= maxDist.num)]
            }
            return(subPairs.gni)
            }) 
    # }
    pairs.gni <- do.call(c,pairs.gni_lst)
    pairs.gni$anchor2.constraint <- NULL
    S4Vectors::mcols(pairs.gni) <- dplyr::rename(data.frame(S4Vectors::mcols(pairs.gni)), constraint="anchor1.constraint")
    names(S4Vectors::mcols(pairs.gni)) <- gsub("anchor1.", "anchor.",gsub("anchor2.", "bait.",names(S4Vectors::mcols(pairs.gni))))
    S4Vectors::mcols(pairs.gni)$name <- paste0(S4Vectors::mcols(pairs.gni)$anchor.bin,"_",S4Vectors::mcols(pairs.gni)$bait.bin)
    S4Vectors::mcols(pairs.gni)$orientation <- (pairs.gni == InteractionSet::swapAnchors(pairs.gni))
    S4Vectors::mcols(pairs.gni)$submatrix.name <- S4Vectors::mcols(pairs.gni)$name
    S4Vectors::mcols(pairs.gni)$submatrix.name[!S4Vectors::mcols(pairs.gni)$orientation] <- paste0(S4Vectors::mcols(pairs.gni)$bait.bin[!S4Vectors::mcols(pairs.gni)$orientation],"_",S4Vectors::mcols(pairs.gni)$anchor.bin[!S4Vectors::mcols(pairs.gni)$orientation])
    columOrder.chr <- unique(c("name","constraint","distance","orientation","submatrix.name","anchor.bin","anchor.name","bait.bin","bait.name",
        names(S4Vectors::mcols(pairs.gni))[grep(x=names(S4Vectors::mcols(pairs.gni)), pattern="^anchor.(.*)$")],
        names(S4Vectors::mcols(pairs.gni))[grep(x=names(S4Vectors::mcols(pairs.gni)), pattern="^bait.(.*)$")]))
    S4Vectors::mcols(pairs.gni) <- dplyr::select(data.frame(S4Vectors::mcols(pairs.gni)) , columOrder.chr) 
    names(pairs.gni) <- S4Vectors::mcols(pairs.gni)$name
    return(pairs.gni)
}