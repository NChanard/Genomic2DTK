#' ExtractSubmatrix
#'
#' Extract matrices in the HiC maps list around genomic features.
#' @param feature.gn <GRanges or Pairs[GRanges] or GInteractions>: The genomic feature on which compute the extraction of HiC submatrix.
#' @param hic.cmx_lst <List[contactMatrix]>: The HiC maps list.
#' @param referencePoint.chr <character>: A character that give the kind of extraction. "rf" to extract regions or "pf" to extract points interactions. (Default "rf")
#' @param res.num <numeric>: the resoulution in used in hic.cmx_lst. If NULL automatically find in resolution attributes of hic.cmx_lst. (Default NULL)
#' @param matriceDim.num <numeric>: The size of matrices in output. (Default 21).
#' @param shiftFactor.num <numeric>: If "referencePoint.chr" is "rf". How much of the distance between anchor and bait is report before and after the region (Default 1).
#' @param cores.num <integer> : An integer to specify the number of cores. (Default 1)
#' @param verbose.bln <logical>: A logical value. If TRUE show the progression in console. (Default TRUE)
#' @return A matrices list.

ExtractSubmatrix <- function(feature.gn=NULL, hic.cmx_lst=NULL, referencePoint.chr="rf",  res.num=NULL, matriceDim.num=21, shiftFactor.num=1,cores.num=1, verbose.bln=TRUE){
        .GInteractionFormatting <- function(feature.gn, res.num){
            if(inherits(feature.gn,"GRanges")){
                feature.gni  <- InteractionSet::GInteractions(GenomicRanges::resize(feature.gn,1,"start"),GenomicRanges::resize(feature.gn,1,"end"))
            }else if(inherits(feature.gn,"Pairs") && inherits(feature.gn@first,"GRanges") && inherits(feature.gn@second,"GRanges")){
                feature.gni <- InteractionSet::GInteractions(feature.gn@first,feature.gn@second)
            }else if (inherits(feature.gn,"GInteractions")){
                feature.gni <- feature.gn
            }
            S4Vectors::mcols(feature.gni) <- S4Vectors::mcols(feature.gn)

            if(is.null(GenomeInfoDb::seqinfo(feature.gni))){
                GenomeInfoDb::seqinfo(feature.gni) <- GenomeInfoDb::seqinfo(feature.gn)
            }
            if(is.null(S4Vectors::mcols(feature.gni)$distance)){
                S4Vectors::mcols(feature.gni)$distance <-  InteractionSet::pairdist(feature.gni)
            }
        
            if(is.null(S4Vectors::mcols(feature.gni)$orientation)){
                S4Vectors::mcols(feature.gni)$orientation <- (feature.gni == InteractionSet::swapAnchors(feature.gni))    
            }
            
            if(is.null(S4Vectors::mcols(feature.gni)$anchor.bin)){
                S4Vectors::mcols(feature.gni)$anchor.bin <-  paste0(GenomeInfoDb::seqnames(InteractionSet::anchors(feature.gni)$first),":", ceiling(InteractionSet::anchors(feature.gni)$first@ranges@start / res.num))

            }    
            
            if(is.null(S4Vectors::mcols(feature.gni)$bait.bin)){
                S4Vectors::mcols(feature.gni)$bait.bin <-  paste0(GenomeInfoDb::seqnames(InteractionSet::anchors(feature.gni)$second),":", ceiling(InteractionSet::anchors(feature.gni)$second@ranges@start / res.num))    
            }

            if(is.null(S4Vectors::mcols(feature.gni)$name)){
                S4Vectors::mcols(feature.gni)$name <- paste0(S4Vectors::mcols(feature.gni)$anchor.bin,"_",S4Vectors::mcols(feature.gni)$bait.bin)
            }
            if(is.null(S4Vectors::mcols(feature.gni)$submatrix.name)){
                S4Vectors::mcols(feature.gni)$submatrix.name <- paste0(S4Vectors::mcols(feature.gni)$anchor.bin,"_",S4Vectors::mcols(feature.gni)$bait.bin)
                S4Vectors::mcols(feature.gni)$submatrix.name[!S4Vectors::mcols(feature.gni)$orientation] <- paste0(S4Vectors::mcols(feature.gni)$bait.bin[!S4Vectors::mcols(feature.gni)$orientation],"_",S4Vectors::mcols(feature.gni)$anchor.bin[!S4Vectors::mcols(feature.gni)$orientation])
            }
            if(!sum(S4Vectors::mcols(feature.gni)$anchor.bin != S4Vectors::mcols(feature.gni)$bait.bin)){
                S4Vectors::mcols(feature.gni)$bait.bin <- NULL
                S4Vectors::mcols(feature.gni)$bin <- S4Vectors::mcols(feature.gni)$anchor.bin
                S4Vectors::mcols(feature.gni)$anchor.bin <- NULL
            }

            return(feature.gni)
        }
    # Run
        # Check Resolution
            if(!is.null(attr(hic.cmx_lst, "resolution"))){res.num <- attr(hic.cmx_lst, "resolution")
            }else if(is.character(res.num)){res.num <- GenomicTK::GenomicSystem(res.num)}
        # Check Dimension
            if(matriceDim.num<5){matriceDim.num<-5}
        # Formatting
            feature.gn <- .GInteractionFormatting(feature.gn=feature.gn, res.num=res.num)
            if(!sum(feature.gn$anchor.bin!=feature.gn$bait.bin)){referencePoint.chr <- "pf"}
        # Resize Features according reference Points
            referencePoint.chr %<>% tolower
            if (referencePoint.chr =="rf"){
                cis.lgk <- DataTK::ReduceRun(
                    GenomeInfoDb::seqnames(InteractionSet::anchors(feature.gn)$first),
                    GenomeInfoDb::seqnames(InteractionSet::anchors(feature.gn)$second),
                    reduceFun.chr="paste",sep="_") %>%
                    as.character %>%
                    lapply(function(combinaison){
                        split <- strsplit(combinaison,"_") %>% unlist 
                        return(split[1]==split[2])
                        }) %>%
                    unlist
                feature.gn <- feature.gn[cis.lgk]
                feature.gn <- feature.gn[which(feature.gn$distance >= (3*res.num))]
                ranges.lst_dtf <- lapply(c("first","second"),function(anchorName.chr){
                    anchor.gnr <- InteractionSet::anchors(feature.gn)[[anchorName.chr]]
                    anchor.dtf <- IRanges::ranges(anchor.gnr) %>%
                        as.data.frame %>%
                        dplyr::select(.data$start,.data$end)
                    magrittr::set_colnames(anchor.dtf,paste0(anchorName.chr,".",names(anchor.dtf))) %>%
                        return(.data)
                    }) 
                ranges.dtf <- do.call(cbind,ranges.lst_dtf) %>%
                    dplyr::mutate(start=pmin(.data$first.start,.data$second.start)) %>%
                    dplyr::mutate(end=pmax(.data$first.end,.data$second.end))
                feature.gnr <- GenomicRanges::GRanges(
                        seqnames = GenomeInfoDb::seqnames(InteractionSet::anchors(feature.gn)$first),
                        ranges = IRanges::IRanges(start=ranges.dtf$start,end=ranges.dtf$end),
                        seqlengths = GenomeInfoDb::seqlengths(feature.gn),
                        seqinfo = GenomeInfoDb::seqinfo(feature.gn)
                    )
                featureResize.gnr <- GenomicRanges::resize(feature.gnr,width=feature.gnr@ranges@width+feature.gn$distance*shiftFactor.num*2,fix="center") 
                featureResize.gni <- InteractionSet::GInteractions(featureResize.gnr, featureResize.gnr)
                S4Vectors::mcols(featureResize.gni ) <- S4Vectors::mcols(feature.gn)
            }else if (referencePoint.chr == "pf") {
                featureResize.gni <- suppressWarnings(GenomicRanges::resize(feature.gn,width=res.num*(matriceDim.num-1)+1,fix="center"))
            }
        # Filt Out Of Bound
            featureFilt.gni <- featureResize.gni[which(
                1L<=data.frame(InteractionSet::anchors(featureResize.gni)$first@ranges)[,"start"] &
                data.frame(InteractionSet::anchors(featureResize.gni)$first@ranges)[,"end"]<=GenomicTK::SeqEnds(InteractionSet::anchors(featureResize.gni)$first) &
                1L<=data.frame(InteractionSet::anchors(featureResize.gni)$second)[,"start"] &
                data.frame(InteractionSet::anchors(featureResize.gni)$second)[,"end"]<=GenomicTK::SeqEnds(InteractionSet::anchors(featureResize.gni)$second)
                )]
        # Filt Duplicated Submatrix before extraction
            featureNoDup.gni <- featureFilt.gni[!duplicated(featureFilt.gni$submatrix.name)]
        # Order according Chromosomes combinaison
            chromosomesCombinaison.rle <- DataTK::ReduceRun(
                GenomeInfoDb::seqnames(InteractionSet::anchors(featureNoDup.gni)$first),
                GenomeInfoDb::seqnames(InteractionSet::anchors(featureNoDup.gni)$second),
                reduceFun.chr="paste",sep="_") 
            order.num <- unlist(lapply(unique(S4Vectors::runValue(chromosomesCombinaison.rle)), function(comb){
                    which(as.character(chromosomesCombinaison.rle)==comb)
            }))
            featureNoDup.gni <- featureNoDup.gni[order.num]
            chromosomesCombinaison.rle <- chromosomesCombinaison.rle[order.num]
        # Separate anchors
            anchors.gnr <- InteractionSet::anchors(featureNoDup.gni)$first
            baits.gnr <- InteractionSet::anchors(featureNoDup.gni)$second
            matAnchors.gnr_lst <- lapply(hic.cmx_lst,InteractionSet::anchors)
        # Extraction
            jobLenght.num <- length(S4Vectors::runValue(chromosomesCombinaison.rle))
            start.tim <- Sys.time()
            submatrix.spm_lst <- do.call(c,lapply(seq_len(jobLenght.num),function(combinaison.ndx){
                combinaisonName.chr <- S4Vectors::runValue(chromosomesCombinaison.rle)[[combinaison.ndx]]
                combinaisonStart.ndx <- cumsum(c(1,S4Vectors::runLength(chromosomesCombinaison.rle)))[[combinaison.ndx]]
                combinaisonEnd.ndx <- cumsum(S4Vectors::runLength(chromosomesCombinaison.rle))[[combinaison.ndx]]
                if(combinaisonName.chr %in% names(matAnchors.gnr_lst)){
                    mat.ndx <- which(names(hic.cmx_lst)==combinaisonName.chr)
                    ovl_row <- GenomicRanges::findOverlaps(anchors.gnr[combinaisonStart.ndx:combinaisonEnd.ndx],matAnchors.gnr_lst[[mat.ndx]]$row) %>%
                        as.data.frame %>%
                        dplyr::group_by(.data$queryHits)
                    ovl_row <- tidyr::nest(ovl_row)
                    ovl_col <- GenomicRanges::findOverlaps(baits.gnr[combinaisonStart.ndx:combinaisonEnd.ndx],matAnchors.gnr_lst[[mat.ndx]]$col) %>%
                        as.data.frame %>%
                        dplyr::group_by(.data$queryHits)
                    ovl_col <- tidyr::nest(ovl_col)
                }else if ({combinaisonName.chr %>% strsplit("_") %>% unlist %>% rev %>% paste(collapse="_")} %in% names(matAnchors.gnr_lst)){
                    mat.ndx <- which(names(hic.cmx_lst)=={combinaisonName.chr %>% strsplit("_") %>% unlist %>% rev %>% paste(collapse="_")})
                    ovl_row <- GenomicRanges::findOverlaps(baits.gnr[combinaisonStart.ndx:combinaisonEnd.ndx],matAnchors.gnr_lst[[mat.ndx]]$row) %>%
                        as.data.frame %>%
                        dplyr::group_by(.data$queryHits)
                    ovl_row <- tidyr::nest(ovl_row)
                    ovl_col <- GenomicRanges::findOverlaps(anchors.gnr[combinaisonStart.ndx:combinaisonEnd.ndx],matAnchors.gnr_lst[[mat.ndx]]$col) %>%
                        as.data.frame %>%
                        dplyr::group_by(.data$queryHits)
                    ovl_col <- tidyr::nest(ovl_col)
                }
                subJobLenght.num <- length(combinaisonStart.ndx:combinaisonEnd.ndx)
                if(cores.num==1){
                    tempSubmatrix.spm_lst <- lapply(seq_len(subJobLenght.num),function(range.ndx){
                        if(verbose.bln){DevTK::ShowLoading(start.tim, range.ndx+(combinaison.ndx-1)*subJobLenght.num,(subJobLenght.num*jobLenght.num))}
                        row.ndx <- unlist(ovl_row[[range.ndx,"data"]],use.names=FALSE)
                        col.ndx <- unlist(ovl_col[[range.ndx,"data"]],use.names=FALSE)
                        if(S4Vectors::metadata(hic.cmx_lst[[mat.ndx]])$type =="cis"){
                            gap.num <- stats::median(col.ndx)-stats::median(row.ndx)
                        }else{
                            gap.num <- Inf
                        }
                        if(gap.num<0){
                            mat.spm <- hic.cmx_lst[[mat.ndx]][col.ndx,row.ndx]@matrix
                        }else{
                            mat.spm <- hic.cmx_lst[[mat.ndx]][row.ndx,col.ndx]@matrix
                        } 
                        if(dim(mat.spm)[1] != matriceDim.num ){
                            mat.spm %<>% DataTK::ResizeMatrix(c(matriceDim.num,matriceDim.num))
                        }
                        if(abs(gap.num) < matriceDim.num ){
                            if(abs(gap.num)>0){
                                mat.spm <- DataTK::PadMtx(mat.mtx=mat.spm, padSize.num=abs(gap.num), value.num=0,side.chr=c("left","bot")) 
                            }
                            mat.spm[lower.tri(mat.spm)] <- NA
                            if(abs(gap.num)>0){
                                mat.spm <- mat.spm[1:matriceDim.num,(abs(gap.num)+1):(matriceDim.num+abs(gap.num))]
                            }
                        }
                        return(as.matrix(mat.spm))
                    })
                }else if(cores.num>=2){
                    parCl <- parallel::makeCluster(cores.num, type ="FORK")
                    doParallel::registerDoParallel(parCl)
                    tempSubmatrix.spm_lst <- parallel::parLapply(parCl,seq_len(subJobLenght.num),function(range.ndx){
                        row.ndx <- unlist(ovl_row[[range.ndx,"data"]],use.names=FALSE)
                        col.ndx <- unlist(ovl_col[[range.ndx,"data"]],use.names=FALSE)
                        if(S4Vectors::metadata(hic.cmx_lst[[mat.ndx]])$type =="cis"){
                            gap.num <- stats::median(col.ndx)-stats::median(row.ndx)
                        }else{
                            gap.num <- Inf
                        }
                        if(gap.num<0){
                            mat.spm <- hic.cmx_lst[[mat.ndx]][col.ndx,row.ndx]@matrix
                        }else{
                            mat.spm <- hic.cmx_lst[[mat.ndx]][row.ndx,col.ndx]@matrix
                        } 
                        if(dim(mat.spm)[1] != matriceDim.num ){
                            mat.spm %<>% DataTK::ResizeMatrix(c(matriceDim.num,matriceDim.num))
                        }
                        if(abs(gap.num) < matriceDim.num ){
                            if(abs(gap.num)>0){
                                mat.spm <- DataTK::PadMtx(mat.mtx=mat.spm, padSize.num=abs(gap.num), value.num=0,side.chr=c("left","bot")) 
                            }
                            mat.spm[lower.tri(mat.spm)] <- NA
                            if(abs(gap.num)>0){
                                mat.spm <- mat.spm[1:matriceDim.num,(abs(gap.num)+1):(matriceDim.num+abs(gap.num))]
                            }
                        }
                        return(as.matrix(mat.spm))
                    })
                    parallel::stopCluster(parCl)
                    DevTK::KillZombies()
                    if(verbose.bln){DevTK::ShowLoading(start.tim,combinaison.ndx, jobLenght.num)}
                }
                return(tempSubmatrix.spm_lst)
            })) %>%
                magrittr::set_names(featureNoDup.gni$submatrix.name) %>%
                magrittr::extract(featureFilt.gni$submatrix.name)
            interactions.ndx <- seq_along(feature.gn$name) %>%
                magrittr::set_names(feature.gn$name) %>%
                magrittr::extract(featureFilt.gni$name)
            attributes(submatrix.spm_lst)$interactions <- feature.gn[interactions.ndx]
            attributes(submatrix.spm_lst)$resolution <- res.num
            attributes(submatrix.spm_lst)$referencePoint <- referencePoint.chr
            attributes(submatrix.spm_lst)$matriceDim <- matriceDim.num
            if(referencePoint.chr=="rf"){attributes(submatrix.spm_lst)$shiftFactor <- shiftFactor.num}
            return(submatrix.spm_lst)
}