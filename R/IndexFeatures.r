#' IndexFeatures
#'
#' Bin multiple GRanges and summary all in one GRange. Could overlap ranges with constraints regions
#' @param gRange.gnr_lst <GRanges or GRangesList or list[GRanges]>: some GRanges or a list of GRanges or a GRangesList
#' @param constraint.gnr <GRanges>: A GRange of constraint regions. If NULL chromosomes in chromSize.dtf are used (Default NULL)
#' @param chromSize.dtf <data.frame>: A data.frame of genome where first colum correspond to the chromosomes names, and the second column correspond to the chromosomes lengths in base pairs.
#' @param binSize.int <integer>: A number that specify the width bins.
#' @param method.chr <character>: A string of a summary method name as 'mean', 'median', 'sum', 'max, 'min'. (Default 'mean'')
#' @param variablesName.chr_vec <character> : A character vector that specify the metadata columns of GRanges on which apply the summary method.
#' @param cores.num <integer> : An integer to specify the number of cores. (Default 1)
#' @param verbose.bln <logical>: A logical value. If TRUE show the progression in console. (Default TRUE)
#' @return A GRange object.
#' @examples
#' library(GenomicED)
#' data("anchors_Peaks.gnr")
#' anchors_Peaks.gnr[1]
#'
#'
#' seqlengths.num <- c('2L'=23513712, '2R'=25286936)
#' chromSize.dtf  <- data.frame(
#'   seqnames   = names(seqlengths.num ), 
#'   seqlengths = seqlengths.num
#'   )
#' binSize.num <- 10000
#'
#' anchors_Index.gnr <- IndexFeatures(
#'   gRange.gnr_lst        = list(Beaf=anchors_Peaks.gnr), 
#'   constraint.gnr        = domains.gnr,
#'   chromSize.dtf         = chromSize.dtf,
#'   binSize.int           = binSize.num,
#'   method.chr            = "max",
#'   variablesName.chr_vec = "score",
#'   cores.num             = 1,
#'   verbose.bln           = FALSE
#'   )
#'
#' anchors_Index.gnr[1]
IndexFeatures <- function(gRange.gnr_lst=NULL, constraint.gnr=NULL, chromSize.dtf=NULL, binSize.int=NULL, method.chr="mean", variablesName.chr_vec=NULL,cores.num=1, verbose.bln=TRUE){
    # Constraint Informations
        if (is.null(constraint.gnr)){
            constraint.gnr <- GenomicRanges::GRanges(
                seqnames = S4Vectors::Rle(chromSize.dtf[, 1]),
                ranges = IRanges::IRanges(start = rep(1, length(chromSize.dtf[, 2])), end = chromSize.dtf[, 2]),
                strand = S4Vectors::Rle('*'),
                name =  chromSize.dtf[, 1]
            )
        }else{
            if(is.null(constraint.gnr$name) | length(which(!is.na(constraint.gnr$name)))==0 ){
                constraint.gnr$name = paste0("Constraint_", seq_along(constraint.gnr))
            }
        }
        seqLevelsStyle.chr <- GenomeInfoDb::seqlevelsStyle(constraint.gnr)
        if(length(seqLevelsStyle.chr)>1){
            seqLevelsStyle.chr <- seqLevelsStyle.chr[[1]]
            GenomeInfoDb::seqlevelsStyle(constraint.gnr) <- seqLevelsStyle.chr
        }
        binnedConstraint.gnr <- GenomicTK::BinGRanges(gRange.gnr=constraint.gnr, chromSize.dtf=chromSize.dtf, binSize.int=binSize.int, verbose.bln=verbose.bln, reduce.bln=FALSE, cores.num=cores.num)
    # Feature Names
        if (inherits(gRange.gnr_lst,"GRanges")){
            gRange.gnr_lst <- list(gRange.gnr_lst) %>% magrittr::set_names("Feature")
        }else if (inherits(gRange.gnr_lst,"GRangesList")){
            gRange.gnr_lst <- as.list(gRange.gnr_lst)
        }
        if (gRange.gnr_lst %>% names %>% is.null){
            gRange.gnr_lst <- magrittr::set_names(gRange.gnr_lst, paste0("Feature_",seq_along(gRange.gnr_lst)))
        }
        gRangeOrder.ndx <- lapply(gRange.gnr_lst,length) %>% unlist %>% order(decreasing=TRUE)
        gRange.gnr_lst <- gRange.gnr_lst[gRangeOrder.ndx]
        feature.chr_vec <- names(gRange.gnr_lst)
    # GRanges Binning
        jobLenght.num <- length(gRange.gnr_lst)
        start.tim <- Sys.time()
        binnedFeature.lst <- lapply(seq_len(jobLenght.num), function(feature.ndx){
            feature.chr <- feature.chr_vec[[feature.ndx]]
            feature.gnr <- IRanges::subsetByOverlaps(gRange.gnr_lst[[feature.chr ]],constraint.gnr)
            GenomeInfoDb::seqlevelsStyle(feature.gnr) <- seqLevelsStyle.chr
            binnedFeature.gnr <- GenomicTK::BinGRanges(gRange.gnr=feature.gnr, chromSize.dtf=chromSize.dtf, binSize.int=binSize.int,  method.chr=method.chr, variablesName.chr_vec=variablesName.chr_vec, verbose.bln=verbose.bln, reduce.bln=TRUE, cores.num=cores.num)
            binnedFeat.tbl <- tibble::tibble(BinnedFeature.ndx = seq_along(binnedFeature.gnr),Feature.name = binnedFeature.gnr$name) %>%
                tidyr::unnest(cols = c(.data$Feature.name)) %>%
                dplyr::group_by(.data$Feature.name)
            binnedFeat.tbl <- tidyr::nest(binnedFeat.tbl) %>%
                magrittr::set_names(c("Feature.name","BinnedFeature.ndx"))
            binnedConstraint.tbl <- tibble::tibble(BinnedConstraint.ndx = seq_along(binnedConstraint.gnr),Constraint.name = binnedConstraint.gnr$name) %>%
                dplyr::group_by(.data$Constraint.name)
            binnedConstraint.tbl <- tidyr::nest(binnedConstraint.tbl) %>%
                magrittr::set_names(c("Constraint.name","BinnedConstraint.ndx"))

            featConstOvlp.ovlp <- GenomicRanges::findOverlaps(feature.gnr, constraint.gnr)
            featConstOvlp.tbl <- tibble::tibble(Feature.name = feature.gnr$name[featConstOvlp.ovlp@from],Constraint.name = constraint.gnr$name[featConstOvlp.ovlp@to]) %>%
                dplyr::left_join(binnedFeat.tbl, by="Feature.name") %>%
                dplyr::select(-(.data$Feature.name)) %>%
                tidyr::unnest(cols = c(.data$BinnedFeature.ndx)) %>%
                unique %>%
                dplyr::group_by(.data$Constraint.name)
            featConstOvlp.tbl <- tidyr::nest(featConstOvlp.tbl) %>%
                magrittr::set_names(c("Constraint.name","BinnedFeature.ndx")) %>%
                dplyr::left_join(binnedConstraint.tbl, by="Constraint.name")
            subJobLenght.num <- featConstOvlp.tbl %>% nrow
            start.tim <- Sys.time()
            if(cores.num==1){
                binnedFeature.gnr_lst <- lapply(seq_len(subJobLenght.num),function(row.ndx){
                    if(verbose.bln){SuperTK::ShowLoading(start.tim, row.ndx+(feature.ndx-1)*subJobLenght.num,(subJobLenght.num*jobLenght.num))}
                    ranges.ndx <- featConstOvlp.tbl$BinnedFeature.ndx[row.ndx] %>% unlist(use.names=FALSE)
                    constraint.ndx <- featConstOvlp.tbl$BinnedConstraint.ndx[row.ndx] %>% unlist(use.names=FALSE)
                    subBinnedFeature.gnr <- IRanges::subsetByOverlaps(binnedFeature.gnr[ranges.ndx], binnedConstraint.gnr[constraint.ndx])
                    subBinnedFeature.gnr$constraint <- featConstOvlp.tbl$Constraint.name[row.ndx]
                    return(subBinnedFeature.gnr)
                })
            }else if(cores.num>=2){
                parCl <- parallel::makeCluster(cores.num, type ="PSOCK")
                parallel::clusterEvalQ(parCl, {
                    library(GenomicRanges)
                })
                binnedFeature.gnr_lst <- parallel::parLapply(parCl,seq_len(subJobLenght.num),function(row.ndx){
                    ranges.ndx <- featConstOvlp.tbl$BinnedFeature.ndx[row.ndx] %>% unlist(use.names=FALSE)
                    constraint.ndx <- featConstOvlp.tbl$BinnedConstraint.ndx[row.ndx] %>% unlist(use.names=FALSE)
                    subBinnedFeature.gnr <- IRanges::subsetByOverlaps(binnedFeature.gnr[ranges.ndx], binnedConstraint.gnr[constraint.ndx])
                    subBinnedFeature.gnr$constraint <- featConstOvlp.tbl$Constraint.name[row.ndx]
                    return(subBinnedFeature.gnr)
                })
                parallel::stopCluster(parCl)
            }
            binnedFeature.gnr <- GenomicTK::MergeGRanges(binnedFeature.gnr_lst, sort.bln=FALSE, reduce.bln=FALSE)
            binnedFeature.gnr$bln <- 1
            names(S4Vectors::mcols(binnedFeature.gnr)) <- paste0(feature.chr, ".",names(S4Vectors::mcols(binnedFeature.gnr)))
            names(S4Vectors::mcols(binnedFeature.gnr))[which(names(S4Vectors::mcols(binnedFeature.gnr)) == paste0(feature.chr, ".bin"))] <- "bin"
            names(S4Vectors::mcols(binnedFeature.gnr))[which(names(S4Vectors::mcols(binnedFeature.gnr)) == paste0(feature.chr, ".constraint"))] <- "constraint"
            binnedFeature.gnr$name <- paste0(binnedFeature.gnr$bin, ":", binnedFeature.gnr$constraint)
            metadataBinnedFeature.dtf <- S4Vectors::mcols(binnedFeature.gnr) %>% data.frame
            S4Vectors::mcols(binnedFeature.gnr) <- NULL
            return(list(binnedFeature.gnr=binnedFeature.gnr, featureMetadata.dtf=metadataBinnedFeature.dtf))
        })

        binnedIndex.gnr <- binnedFeature.lst %>% lapply("[[", "binnedFeature.gnr") %>% GenomicTK::MergeGRanges(sort.bln=FALSE, reduce.bln=FALSE)
        featureMetadata.lst_dtf <- binnedFeature.lst %>% lapply("[[", "featureMetadata.dtf")
        S4Vectors::mcols(binnedIndex.gnr) <- SuperTK::BindFillRows(featureMetadata.lst_dtf)
        ids.lst <- binnedIndex.gnr$name
        dupplicatedIds.lst <- unique(ids.lst[duplicated(ids.lst)])
        idDuplicated.ndx <- which(ids.lst %in% dupplicatedIds.lst)
    # Merge GRanges into one index and duplicated Bin handle
        if(length(idDuplicated.ndx)){
            binnedIndexDuplicated.dtf <- data.frame(binnedIndex.gnr[idDuplicated.ndx])
            binnedIndexDuplicated.tbl <- tibble::tibble(binnedIndexDuplicated.dtf) %>% dplyr::group_by(.data$name) 
            binnedIndexDuplicated.tbl <- tidyr::nest(binnedIndexDuplicated.tbl)
            binnedIndexNoDuplicated.dtf <- data.frame(binnedIndex.gnr[-idDuplicated.ndx])
            binnedIndexNoDuplicated.tbl <- tibble::tibble(binnedIndexNoDuplicated.dtf)
            start.tim <- Sys.time()
            jobLenght.num <- nrow(binnedIndexDuplicated.tbl)
            if(cores.num==1){
                binnedIndexDuplicated.lst <- lapply(seq_len(jobLenght.num), function(row.ndx){
                    if(verbose.bln){SuperTK::ShowLoading(start.tim, row.ndx, jobLenght.num)}
                    rowName.chr <- binnedIndexDuplicated.tbl$name[[row.ndx]]
                    row <- binnedIndexDuplicated.tbl$data[[row.ndx]]
                    col.lst <- lapply(seq_along(row),function(col.ndx){
                        col <- dplyr::pull(row,col.ndx)
                        if(length(unique(stats::na.omit(col)))==0){
                            return(NA)
                        }else if(length(unique(stats::na.omit(col)))==1) {
                            return(unique(stats::na.omit(col)))
                        }else {
                            return(list(unlist(col)))
                            }
                        })  %>% magrittr::set_names(names(row))

                    tibble::as_tibble(col.lst) %>%
                    tibble::add_column(name = rowName.chr) %>%
                    return(.data)
                })
            }else if(cores.num>=2){
                parCl <- parallel::makeCluster(cores.num, type ="PSOCK")
                binnedIndexDuplicated.lst <- parallel::parLapply(parCl,seq_len(jobLenght.num), function(row.ndx){
                    rowName.chr <- binnedIndexDuplicated.tbl$name[[row.ndx]]
                    row <- binnedIndexDuplicated.tbl$data[[row.ndx]]
                    col.lst <- lapply(seq_along(row),function(col.ndx){
                        col <- dplyr::pull(row,col.ndx)
                        if(length(unique(stats::na.omit(col)))==0){
                            return(NA)
                        }else if(length(unique(stats::na.omit(col)))==1) {
                            return(unique(stats::na.omit(col)))
                        }else {
                            return(list(unlist(col)))
                        }
                        })  %>% magrittr::set_names(names(row))
                    tibble::as_tibble(col.lst) %>%
                    tibble::add_column(name = rowName.chr) %>%
                    return(.data)
                })
                parallel::stopCluster(parCl)
            }
            binnedIndexDuplicated.tbl <- dplyr::bind_rows(binnedIndexDuplicated.lst)
            binnedIndex.gnr <- rbind(binnedIndexDuplicated.tbl, binnedIndexNoDuplicated.tbl) %>% data.frame %>% methods::as("GRanges")
        }
        for(featureName.chr in feature.chr_vec){
            colname.chr =  paste0(featureName.chr, ".bln")
            S4Vectors::mcols(binnedIndex.gnr)[which(is.na(S4Vectors::mcols(binnedIndex.gnr)[, colname.chr])),colname.chr] <- 0
            S4Vectors::mcols(binnedIndex.gnr)[,colname.chr] <- methods::as(S4Vectors::mcols(binnedIndex.gnr)[, colname.chr], "Rle")
        }
        S4Vectors::mcols(binnedIndex.gnr) %<>% as.data.frame %>% dplyr::select(.data$name,.data$bin,.data$constraint,tidyselect::everything()) 
    return(sort(binnedIndex.gnr))
}