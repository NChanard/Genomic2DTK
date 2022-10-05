#' GnpToCm
#'
#' Convert a Pairs of GRanges in an list of ContactMatrix (matricies are stored in sparseMatrix format).
#' @param hic.gnp <Pairs[GRanges]>: some GRanges or a list of GRanges or a GRangesList
#' @param res.num <integer>: A number resolution of the HiC.
#' @param chromSize.dtf <data.frame>: A data.frame of genome where first colum correspond to the chromosomes names, and the second column correspond to the chromosomes lengths in base pairs.
#' @param cores.num <integer> : An integer to specify the number of cores. (Default 1)
#' @param verbose.bln <logical>: A logical value. If TRUE show the progression in console. (Default TRUE)
#' @return A list of ContactMatrix.

GnpToCm <- function(hic.gnp=NULL, res.num=NULL, chromSize.dtf=NULL, verbose.bln=TRUE, cores.num=1){
        seqlengths.lst <- chromSize.dtf %>%
            dplyr::pull(2) %>%
            magrittr::set_names({chromSize.dtf%>% dplyr::pull(1)})
        binnedGenome.grn <- GenomicRanges::tileGenome(seqlengths.lst, tilewidth=res.num, cut.last.tile.in.chrom=TRUE)
    # Chromosomes Size
        chromSize.tbl <- tibble::tibble(magrittr::set_colnames(chromSize.dtf,c("seqnames","seqlengths")))
        chromSize.tbl$dimension <- chromSize.tbl %>% magrittr::extract("seqlengths") %>% magrittr::divide_by(res.num) %>% ceiling %>% unlist %>% as.numeric
        genome.gnr = GenomicRanges::GRanges(
                seqnames = S4Vectors::Rle(chromSize.dtf[, 1]),
                ranges = IRanges::IRanges(start = rep(1, length(chromSize.dtf[, 2])), end = chromSize.dtf[, 2]),
                strand = S4Vectors::Rle('*'),
                name =  S4Vectors::Rle(chromSize.dtf[, 1])
            )
        names(genome.gnr) <- genome.gnr$name
    # Combinaison present in bedpe
        chromosomesCombinaison.rle = SuperTK::ReduceRun(GenomicRanges::seqnames(hic.gnp@first),GenomicRanges::seqnames(hic.gnp@second),reduce.fun="paste",sep="_")
        chromosomesCombinaison.chr = S4Vectors::runValue(chromosomesCombinaison.rle)
    # Combinaison present in bedpe
        chromosomesCombinaison.chr <- unique(paste0(hic.gnp@first@seqnames, "_", hic.gnp@second@seqnames))
        matrixSymmetric.bln <- chromosomesCombinaison.chr %>%
            strsplit("_") %>%
            lapply(function(name.chr){name.chr[[1]]  == name.chr[[2]]}) %>%
            unlist
        matrixType.str <- chromosomesCombinaison.chr
        matrixType.str[which(matrixSymmetric.bln)] <- "cis"
        matrixType.str[which(!matrixSymmetric.bln)] <- "trans"
        matrixKind.str <- chromosomesCombinaison.chr
        matrixKind.str[which(matrixSymmetric.bln)] <- "U"
        matrixKind.str[which(!matrixSymmetric.bln)] <- NA
        attributes.tbl <-  dplyr::bind_cols(name=chromosomesCombinaison.chr, type=matrixType.str, kind=matrixKind.str, symmetric=matrixSymmetric.bln )
    # Convert
        jobLenght.num <- length(chromosomesCombinaison.chr)
        if(cores.num==1){
            start.tim <- Sys.time()
            if(verbose.bln){cat("\n")}
            hic.cmx_lst <- lapply(seq_len(jobLenght.num),function(combin.ndx){
                if(verbose.bln){SuperTK::ShowLoading(start.tim,combin.ndx,jobLenght.num)}
                combin.lst <- chromosomesCombinaison.chr[[combin.ndx]] %>%
                    strsplit("_") %>%
                    unlist
                start.ndx <- c(1,S4Vectors::runLength(chromosomesCombinaison.rle)) %>%
                    cumsum %>%
                    magrittr::extract(combin.ndx)
                end.ndx <- S4Vectors::runLength(chromosomesCombinaison.rle) %>%
                    cumsum %>%
                    magrittr::extract(combin.ndx)
                i.num_vec <- hic.gnp@first@ranges@start %>%
                    magrittr::extract(start.ndx:end.ndx) %>%
                    magrittr::divide_by(res.num) %>%
                    ceiling %>%
                    as.integer
                j.num_vec <- hic.gnp@second@ranges@start %>%
                    magrittr::extract(start.ndx:end.ndx) %>%
                    magrittr::divide_by(res.num) %>%
                    ceiling %>%
                    as.integer
                x.num_vec <- S4Vectors::mcols(hic.gnp) %>%
                    magrittr::use_series("score") %>%
                    magrittr::extract(start.ndx:end.ndx)
                dimensionI.ndx <- chromSize.tbl %>%
                    dplyr::pull("seqnames") %>%
                    magrittr::equals(combin.lst[1]) %>%
                    which
                dimensionI.num <- (chromSize.tbl %>%  dplyr::pull("dimension"))[dimensionI.ndx]
                dimensionJ.ndx <- chromSize.tbl %>%
                    dplyr::pull("seqnames") %>%
                    magrittr::equals(combin.lst[2]) %>%
                    which
                dimensionJ.num <- (chromSize.tbl %>%  dplyr::pull("dimension"))[dimensionJ.ndx]
                attributes.lst <- attributes.tbl[combin.ndx,] %>%
                    tibble::add_column(resolution = res.num) %>%
                    as.list
                hic.spm <- Matrix::sparseMatrix(i=i.num_vec, j=j.num_vec, x=x.num_vec, dims=c(dimensionI.num,dimensionJ.num))
                row.regions = binnedGenome.grn[GenomicRanges::findOverlaps(binnedGenome.grn,genome.gnr[combin.lst[[1]]])@from]
                col.regions = binnedGenome.grn[GenomicRanges::findOverlaps(binnedGenome.grn,genome.gnr[combin.lst[[2]]])@from]
                hic.cmx = InteractionSet::ContactMatrix(hic.spm, row.regions, col.regions) 
                S4Vectors::metadata(hic.cmx) <- attributes.lst
                return(hic.cmx)
            })
            if(verbose.bln){cat("\n")}
        }else if(cores.num>=2){
            parCl <- parallel::makeCluster(cores.num, type ="FORK")
            doParallel::registerDoParallel(parCl)
            hic.cmx_lst <- parallel::parLapply(parCl,seq_len(jobLenght.num),function(combin.ndx){
                combin.lst <- chromosomesCombinaison.chr[[combin.ndx]] %>%
                    strsplit("_") %>%
                    unlist
                start.ndx <- c(1,S4Vectors::runLength(chromosomesCombinaison.rle)) %>%
                    cumsum %>%
                    magrittr::extract(combin.ndx)
                end.ndx <- S4Vectors::runLength(chromosomesCombinaison.rle) %>%
                    cumsum %>%
                    magrittr::extract(combin.ndx)
                i.num_vec <- hic.gnp@first@ranges@start %>%
                    magrittr::extract(start.ndx:end.ndx) %>%
                    magrittr::divide_by(res.num) %>%
                    ceiling %>%
                    as.integer
                j.num_vec <- hic.gnp@second@ranges@start %>%
                    magrittr::extract(start.ndx:end.ndx) %>%
                    magrittr::divide_by(res.num) %>%
                    ceiling %>%
                    as.integer
                x.num_vec <- S4Vectors::mcols(hic.gnp) %>%
                    magrittr::use_series("score") %>%
                    magrittr::extract(start.ndx:end.ndx)
                dimensionI.ndx <- chromSize.tbl %>%
                    dplyr::pull("seqnames") %>%
                    magrittr::equals(combin.lst[1]) %>%
                    which
                dimensionI.num <- (chromSize.tbl %>%  dplyr::pull("dimension"))[dimensionI.ndx]
                dimensionJ.ndx <- chromSize.tbl %>%
                    dplyr::pull("seqnames") %>%
                    magrittr::equals(combin.lst[2]) %>%
                    which
                dimensionJ.num <- (chromSize.tbl %>%  dplyr::pull("dimension"))[dimensionJ.ndx]
                attributes.lst <- attributes.tbl[combin.ndx,] %>%
                    tibble::add_column(resolution = res.num) %>%
                    as.list
                hic.spm <- Matrix::sparseMatrix(i=i.num_vec, j=j.num_vec, x=x.num_vec, dims=c(dimensionI.num,dimensionJ.num))
                row.regions = binnedGenome.grn[GenomicRanges::findOverlaps(binnedGenome.grn,genome.gnr[combin.lst[[1]]])@from]
                col.regions = binnedGenome.grn[GenomicRanges::findOverlaps(binnedGenome.grn,genome.gnr[combin.lst[[2]]])@from]
                hic.cmx = InteractionSet::ContactMatrix(hic.spm, row.regions, col.regions) 
                S4Vectors::metadata(hic.cmx) <- attributes.lst
                return(hic.cmx)
            })
            parallel::stopCluster(parCl)
        }
        hic.cmx_lst %>%
        magrittr::set_names(chromosomesCombinaison.chr) %>%
        SuperTK::AddAttr(list(resolution=res.num, chromSize = chromSize.tbl, matricesKind=attributes.tbl )) %>%
        return(.data)
}