#' Function that indexes a GRanges object on binned genome and constraints. Needed prior SearchPairs() function. Note that bins in the same constraint region only will be paired in SearchPairs().
#' 
#' IndexFeatures
#' @description Function that indexes a GRanges object on binned genome and constraints. Needed prior SearchPairs() function.
#' @param gRange.gnr_lst <GRanges or GRangesList or list[GRanges]>: GRanges object, list of GRanges or GRangesList containing coordinates to index.
#' @param constraint.gnr <GRanges>: GRanges object of constraint regions. Note that bins in the same constraint region only will be paired in SearchPairs(). If NULL chromosomes in chromSize.dtf are used as constraints (Default NULL)
#' @param chromSize.dtf <data.frame>: A data.frame containing chromosomes names and lengths in base pairs (see example).
#' @param binSize.num <integer>: Bin size in bp - corresponds to HiC matrix resolution.
#' @param variablesName.chr_vec <character> : A character vector that specify the metadata columns of GRanges on which apply the summary method if multiple ranges are indexed in the same bin.
#' @param method.chr <character>: A string defining which summary method is used on metadata columns defined in variablesName.chr_vec if multiple ranges are indexed in the same bin. Use 'mean', 'median', 'sum', 'max' or 'min'. (Default 'mean'')
#' @param cores.num <integer> : Number of cores used. (Default 1)
#' @param verbose.bln <logical>: A logical value. If TRUE show the progression in console. (Default TRUE)
#' @return A GRanges object.
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
#'   constraint.gnr        = NULL,
#'   chromSize.dtf         = chromSize.dtf,
#'   binSize.num           = binSize.num,
#'   method.chr            = "max",
#'   variablesName.chr_vec = "score",
#'   cores.num             = 1,
#'   verbose.bln           = FALSE
#'   )
#'
#' anchors_Index.gnr[1]
IndexFeatures <- function(gRange.gnr_lst=NULL, constraint.gnr=NULL, chromSize.dtf=NULL, binSize.num=NULL, method.chr="mean", variablesName.chr_vec=NULL,cores.num=1, verbose.bln=TRUE){
    # Constraint Informations
        if (is.null(constraint.gnr)){
            constraint.gnr <- GenomicRanges::GRanges(
                seqnames = chromSize.dtf[, 1],
                ranges = IRanges::IRanges(start = rep(1, length(chromSize.dtf[, 2])), end = chromSize.dtf[, 2]),
                strand = '*',
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
        binnedConstraint.gnr <- GenomicTK::BinGRanges(gRange.gnr=constraint.gnr, chromSize.dtf=chromSize.dtf, binSize.num=binSize.num, verbose.bln=verbose.bln, reduce.bln=FALSE, cores.num=cores.num)
    # Feature Names
        if (inherits(gRange.gnr_lst,"GRanges")){
            gRange.gnr_lst <- list(Features = gRange.gnr_lst)
        }else if (inherits(gRange.gnr_lst,"GRangesList")){
            gRange.gnr_lst <- as.list(gRange.gnr_lst)
        }
        if (is.null(names(gRange.gnr_lst))){
            gRange.gnr_lst <- stats::setNames(gRange.gnr_lst, paste0("Feature_",seq_along(gRange.gnr_lst)))
        }
        gRangeOrder.ndx <- lapply(gRange.gnr_lst,length) |> unlist() |> order(decreasing=TRUE)
        gRange.gnr_lst <- gRange.gnr_lst[gRangeOrder.ndx]
        feature.chr_vec <- names(gRange.gnr_lst)
    # GRanges Binning
        jobLenght.num <- length(gRange.gnr_lst)
        start.tim <- Sys.time()
        binnedFeature.lst <- lapply(seq_len(jobLenght.num), function(feature.ndx){
            feature.chr <- feature.chr_vec[[feature.ndx]]
            feature.gnr <- IRanges::subsetByOverlaps(gRange.gnr_lst[[feature.chr ]],constraint.gnr)
            GenomeInfoDb::seqlevelsStyle(feature.gnr) <- seqLevelsStyle.chr
            binnedFeature.gnr <- GenomicTK::BinGRanges(gRange.gnr=feature.gnr, chromSize.dtf=chromSize.dtf, binSize.num=binSize.num,  method.chr=method.chr, variablesName.chr_vec=variablesName.chr_vec, verbose.bln=verbose.bln, reduce.bln=TRUE, cores.num=cores.num)
            binnedFeat.tbl <- tibble::tibble(BinnedFeature.ndx = seq_along(binnedFeature.gnr),Feature.name = binnedFeature.gnr$name) |>
                tidyr::unnest(cols = "Feature.name")
            binnedFeat.tbl <- dplyr::group_by(binnedFeat.tbl, Feature.name = binnedFeat.tbl$Feature.name)
            binnedFeat.tbl <- tidyr::nest(binnedFeat.tbl) |>
                stats::setNames(c("Feature.name","BinnedFeature.ndx"))
            binnedConstraint.tbl <- tibble::tibble(BinnedConstraint.ndx = seq_along(binnedConstraint.gnr),Constraint.name = binnedConstraint.gnr$name) 
            binnedConstraint.tbl <- dplyr::group_by(binnedConstraint.tbl, Constraint.name = binnedConstraint.tbl$Constraint.name)
            binnedConstraint.tbl <- tidyr::nest(binnedConstraint.tbl) |>
                stats::setNames(c("Constraint.name","BinnedConstraint.ndx"))
            featConstOvlp.ovlp <- GenomicRanges::findOverlaps(feature.gnr, constraint.gnr)
            featConstOvlp.tbl <- tibble::tibble(Feature.name = feature.gnr$name[featConstOvlp.ovlp@from],Constraint.name = constraint.gnr$name[featConstOvlp.ovlp@to]) |>
                dplyr::left_join(binnedFeat.tbl, by="Feature.name") |>
                dplyr::select(-"Feature.name") |>
                tidyr::unnest(cols = "BinnedFeature.ndx") |>
                unique()
            featConstOvlp.tbl <- dplyr::group_by(featConstOvlp.tbl, Constraint.name = featConstOvlp.tbl$Constraint.name)
            featConstOvlp.tbl <- tidyr::nest(featConstOvlp.tbl) |>
                stats::setNames(c("Constraint.name","BinnedFeature.ndx")) |>
                dplyr::left_join(binnedConstraint.tbl, by="Constraint.name")
            subJobLenght.num <- nrow(featConstOvlp.tbl)
            start.tim <- Sys.time()
            if(cores.num==1){
                binnedFeature.gnr_lst <- lapply(seq_len(subJobLenght.num),function(row.ndx){
                    if(verbose.bln){SuperTK::ShowLoading(start.tim, row.ndx+(feature.ndx-1)*subJobLenght.num,(subJobLenght.num*jobLenght.num))}
                    ranges.ndx <- featConstOvlp.tbl$BinnedFeature.ndx[row.ndx] |> unlist(use.names=FALSE)
                    constraint.ndx <- featConstOvlp.tbl$BinnedConstraint.ndx[row.ndx] |> unlist(use.names=FALSE)
                    subBinnedFeature.gnr <- IRanges::subsetByOverlaps(binnedFeature.gnr[ranges.ndx], binnedConstraint.gnr[constraint.ndx])
                    subBinnedFeature.gnr$constraint <- featConstOvlp.tbl$Constraint.name[row.ndx]
                    return(subBinnedFeature.gnr)
                })
            }else if(cores.num>=2){
                multicoreParam <- BiocParallel::MulticoreParam(workers = cores.num) # DD221108 change to BiocParallel
                # parCl <- parallel::makeCluster(cores.num, type ="PSOCK")
                # parallel::clusterEvalQ(parCl, {
                #     library(GenomicRanges)
                # })
                # binnedFeature.gnr_lst <- parallel::parLapply(parCl,seq_len(subJobLenght.num),function(row.ndx){
                binnedFeature.gnr_lst <- BiocParallel::bplapply(BPPARAM = multicoreParam,seq_len(subJobLenght.num),function(row.ndx){
                    ranges.ndx <- featConstOvlp.tbl$BinnedFeature.ndx[row.ndx] |> unlist(use.names=FALSE)
                    constraint.ndx <- featConstOvlp.tbl$BinnedConstraint.ndx[row.ndx] |> unlist(use.names=FALSE)
                    subBinnedFeature.gnr <- IRanges::subsetByOverlaps(binnedFeature.gnr[ranges.ndx], binnedConstraint.gnr[constraint.ndx])
                    subBinnedFeature.gnr$constraint <- featConstOvlp.tbl$Constraint.name[row.ndx]
                    return(subBinnedFeature.gnr)
                })
                # parallel::stopCluster(parCl)
            }
            binnedFeature.gnr <- GenomicTK::MergeGRanges(binnedFeature.gnr_lst, sort.bln=FALSE, reduce.bln=FALSE)
            # binnedFeature.gnr$bln <- 1
            binnedFeature.gnr$bln <- T # DD221108 change for T/F instead of 1/0
            names(S4Vectors::mcols(binnedFeature.gnr)) <- paste0(feature.chr, ".",names(S4Vectors::mcols(binnedFeature.gnr)))
            names(S4Vectors::mcols(binnedFeature.gnr))[which(names(S4Vectors::mcols(binnedFeature.gnr)) == paste0(feature.chr, ".bin"))] <- "bin"
            names(S4Vectors::mcols(binnedFeature.gnr))[which(names(S4Vectors::mcols(binnedFeature.gnr)) == paste0(feature.chr, ".constraint"))] <- "constraint"
            binnedFeature.gnr$name <- paste0(binnedFeature.gnr$bin, ":", binnedFeature.gnr$constraint)
            metadataBinnedFeature.dtf <- data.frame(S4Vectors::mcols(binnedFeature.gnr))
            S4Vectors::mcols(binnedFeature.gnr) <- NULL
            return(list(binnedFeature.gnr=binnedFeature.gnr, featureMetadata.dtf=metadataBinnedFeature.dtf))
        })

        binnedIndex.gnr <- binnedFeature.lst |> lapply("[[", "binnedFeature.gnr") |> GenomicTK::MergeGRanges(sort.bln=FALSE, reduce.bln=FALSE)
        featureMetadata.lst_dtf <- binnedFeature.lst |> lapply("[[", "featureMetadata.dtf")
        S4Vectors::mcols(binnedIndex.gnr) <- SuperTK::BindFillRows(featureMetadata.lst_dtf)
        ids.lst <- binnedIndex.gnr$name
        dupplicatedIds.lst <- unique(ids.lst[duplicated(ids.lst)])
        idDuplicated.ndx <- which(ids.lst %in% dupplicatedIds.lst)
    # Merge GRanges into one index and duplicated Bin handle
        if(length(idDuplicated.ndx)){
            binnedIndexDuplicated.dtf <- data.frame(binnedIndex.gnr[idDuplicated.ndx])
            binnedIndexDuplicated.tbl <- tibble::tibble(binnedIndexDuplicated.dtf) 
            binnedIndexDuplicated.tbl <- dplyr::group_by(binnedIndexDuplicated.tbl, name = binnedIndexDuplicated.tbl$name)
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
                        })  |> stats::setNames(names(row))
                    binnedIndexDuplicated.tbl <- tibble::as_tibble(col.lst) |>
                        tibble::add_column(name = rowName.chr)
                    return(binnedIndexDuplicated.tbl)
                })
            }else if(cores.num>=2){
                multicoreParam <- BiocParallel::MulticoreParam(workers = cores.num) # DD221108 change to BiocParallel
                binnedIndexDuplicated.lst <- BiocParallel::bplapply(BPPARAM = multicoreParam,seq_len(jobLenght.num), function(row.ndx){
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
                        })  |> stats::setNames(names(row))
                    binnedIndexDuplicated.tbl <- tibble::as_tibble(col.lst) |>
                        tibble::add_column(name = rowName.chr)
                    return(binnedIndexDuplicated.tbl)
                })
                # parallel::stopCluster(parCl)
            }
            binnedIndexDuplicated.tbl <- dplyr::bind_rows(binnedIndexDuplicated.lst)
            binnedIndex.gnr <- rbind(binnedIndexDuplicated.tbl, binnedIndexNoDuplicated.tbl) |> data.frame() |> methods::as("GRanges")
        }
        for(featureName.chr in feature.chr_vec){
            colname.chr =  paste0(featureName.chr, ".bln")
            # S4Vectors::mcols(binnedIndex.gnr)[which(is.na(S4Vectors::mcols(binnedIndex.gnr)[, colname.chr])),colname.chr] <- 0
            S4Vectors::mcols(binnedIndex.gnr)[which(is.na(S4Vectors::mcols(binnedIndex.gnr)[, colname.chr])),colname.chr] <- F # DD221108 change for T/F instead of 1/0
            S4Vectors::mcols(binnedIndex.gnr)[,colname.chr] <- methods::as(S4Vectors::mcols(binnedIndex.gnr)[, colname.chr], "Rle")
        }
        columOrder.chr <- unique(c("name","bin","constraint",names(S4Vectors::mcols(binnedIndex.gnr))))
        S4Vectors::mcols(binnedIndex.gnr) <- as.data.frame(S4Vectors::mcols(binnedIndex.gnr)) |> dplyr::select(columOrder.chr) 
    return(sort(binnedIndex.gnr))
}