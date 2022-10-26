#' Cut HiC map in HiC chunks.
#' 
#' CutHiC
#' @description Cut a mega contactMatrix (join of multiple chromosomic maps) inq a list of contactMatrix.
#' @param megaHic.cmx <contactMatrix>: The HiC megamap.
#' @param cores.num <numerical> : An integer to specify the number of cores. (Default 1)
#' @param verbose.bln <logical>: A logical value. If TRUE show the progression in console. (Default TRUE)
#' @return A matrices list.
#' @examples
CutHiC <- function(megaHic.cmx, cores.num=1, verbose.bln=TRUE){
    res.num <- megaHic.cmx@metadata$resolution
    mtx.chr <- megaHic.cmx@metadata$mtx
    chromSize.dtf <-  megaHic.cmx@metadata$chromSize
    binnedGenome.grn <- chromSize.dtf |>
        dplyr::pull("length") |>
        stats::setNames(chromSize.dtf$name) |>
        GenomicRanges::tileGenome(tilewidth=res.num, cut.last.tile.in.chrom=TRUE)
    GenomeInfoDb::seqlengths(binnedGenome.grn) <- chromSize.dtf$length |> stats::setNames(chromSize.dtf$name) 
    attributes.tbl <- megaHic.cmx@metadata$matricesKind
    chromComb.lst <- attributes.tbl$name
    if(cores.num==1){
        start.tim <- Sys.time()
        jobLenght.num <- length(chromComb.lst)

        hic.lst_cmx <- lapply(seq_along(chromComb.lst), function(ele.ndx){
            if(verbose.bln){SuperTK::ShowLoading(start.tim, ele.ndx,jobLenght.num)}
            # Chromosomes
                ele.lst <- unlist(strsplit(chromComb.lst[[ele.ndx]],"_"))
                row.regions = binnedGenome.grn[which(as.vector(binnedGenome.grn@seqnames) ==  ele.lst[[1]])]
                col.regions = binnedGenome.grn[which(as.vector(binnedGenome.grn@seqnames) == ele.lst[[2]])]
            # Extract matrix
                hic.spm <- megaHic.cmx@matrix[GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$row, row.regions)@from,GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$column, col.regions)@from]
                hic.cmx <- InteractionSet::ContactMatrix(hic.spm, row.regions, col.regions)
            # removedCounts
                removedCounts <- NULL
                if(!is.null(megaHic.cmx@metadata$removedCounts)){
                    removedCounts <- list(removedCounts=megaHic.cmx@metadata$removedCounts[GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$row, row.regions)@from,GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$column, col.regions)@from])
                }
            # observed
                observed <- NULL
                if(!is.null(megaHic.cmx@metadata$observed)){
                    observed.spm <- megaHic.cmx@matrix
                    observed.spm@x <- megaHic.cmx@metadata$observed
                    observed <- list(observed=observed.spm[GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$row, row.regions)@from,GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$column, col.regions)@from]@x)
                }
            # normalizer
                normalizer <- NULL
                if(!is.null(megaHic.cmx@metadata$normalizer)){
                    normalizer.spm <- megaHic.cmx@matrix
                    normalizer.spm@x <- megaHic.cmx@metadata$normalizer  
                    normalizer <- list(normalizer=normalizer.spm[GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$row, row.regions)@from,GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$column, col.regions)@from]@x)
                }
            # expected
                expected <- NULL
                if(!is.null(megaHic.cmx@metadata$expected)){
                    expected.spm <- megaHic.cmx@matrix
                    expected.spm@x <- megaHic.cmx@metadata$expected  
                    expected <- list(expected=expected.spm[GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$row, row.regions)@from,GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$column, col.regions)@from]@x)
                }
            # Metadata
                hic.cmx@metadata <- dplyr::filter(attributes.tbl, attributes.tbl$name == chromComb.lst[[ele.ndx]]) |>
                        tibble::add_column(resolution = res.num) |>
                        as.list() |> c(removedCounts,observed,normalizer,expected)
            return(hic.cmx)
        })
    } else if(cores.num>=2){
        parCl <- parallel::makeCluster(cores.num, type ="FORK")
            hic.lst_cmx <- parallel::parLapply(parCl,seq_along(chromComb.lst), function(ele.ndx){
                # Chromosomes
                    ele.lst <- unlist(strsplit(chromComb.lst[[ele.ndx]],"_"))
                    row.regions = binnedGenome.grn[which(as.vector(binnedGenome.grn@seqnames) ==  ele.lst[[1]])]
                    col.regions = binnedGenome.grn[which(as.vector(binnedGenome.grn@seqnames) == ele.lst[[2]])]
                # Extract matrix
                    hic.spm <- megaHic.cmx@matrix[GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$row, row.regions)@from,GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$column, col.regions)@from]
                    hic.cmx <- InteractionSet::ContactMatrix(hic.spm, row.regions, col.regions)
                # removedCounts
                    removedCounts <- NULL
                    if(!is.null(megaHic.cmx@metadata$removedCounts)){
                        removedCounts <- list(removedCounts=megaHic.cmx@metadata$removedCounts[GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$row, row.regions)@from,GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$column, col.regions)@from])
                    }
                # observed
                    observed <- NULL
                    if(!is.null(megaHic.cmx@metadata$observed)){
                        observed.spm <- megaHic.cmx@matrix
                        observed.spm@x <- megaHic.cmx@metadata$observed
                        observed <- list(observed=observed.spm[GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$row, row.regions)@from,GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$column, col.regions)@from]@x)
                    }
                # normalizer
                    normalizer <- NULL
                    if(!is.null(megaHic.cmx@metadata$normalizer)){
                        normalizer.spm <- megaHic.cmx@matrix
                        normalizer.spm@x <- megaHic.cmx@metadata$normalizer  
                        normalizer <- list(normalizer=normalizer.spm[GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$row, row.regions)@from,GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$column, col.regions)@from]@x)
                    }
                # expected
                    expected <- NULL
                    if(!is.null(megaHic.cmx@metadata$expected)){
                        expected.spm <- megaHic.cmx@matrix
                        expected.spm@x <- megaHic.cmx@metadata$expected  
                        expected <- list(expected=expected.spm[GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$row, row.regions)@from,GenomicRanges::findOverlaps(InteractionSet::anchors(megaHic.cmx)$column, col.regions)@from]@x)
                    }
                # Metadata
                    hic.cmx@metadata <- dplyr::filter(attributes.tbl, attributes.tbl$name == chromComb.lst[[ele.ndx]]) |>
                            tibble::add_column(resolution = res.num) |>
                            as.list() |> c(removedCounts,observed,normalizer)
                return(hic.cmx)
            })
        parallel::stopCluster(parCl)
    }
    hic.lst_cmx <- hic.lst_cmx |>
        stats::setNames(chromComb.lst) |>
        SuperTK::AddAttr(list(
            resolution = res.num,
            mtx = mtx.chr,
            chromSize = tibble::as_tibble(chromSize.dtf),
            matricesKind=attributes.tbl
            ))

    return(hic.lst_cmx)
}