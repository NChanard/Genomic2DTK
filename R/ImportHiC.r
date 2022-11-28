#' Import Hic data 
#'
#' ImportHiC
#' @description Import ..hic, .cool, .mcool or .bedpe data  
#' @param file.pth <GRanges or Pairs[GRanges] or GInteractions>: The genomic feature on which compute the extraction of HiC submatrix. Extension should be .hic, .cool, .mcool, .h5, .hdf5, .HDF5 or .bedpe" assuming .h5 et .hdf5 are only for cool (not mcool).
#' @param res.num <numeric>: the HiC resolution.
#' @param chromSize.dtf <data.frame>: a data.frame where first colum correspond to the chromosomes names, and the second column correspond to the chromosomes lengths in base pairs.
#' @param chrom_1.chr <numeric>: The seqnames of firsts chromosmes (rows in matrix).
#' @param chrom_2.chr <numeric>: The seqnames of second chromosmes (col in matrix). If is NULL is equal to chrom_1.chr (Defalt NULL)
#' @param cores.num <numerical> : An integer to specify the number of cores. (Default 1)
#' @param verbose.bln <logical>: A logical value. If TRUE show the progression in console. (Default TRUE)
#' @return A matrices list.
#' @examples
#' \dontrun{
#' # Download HiC File
#' options(timeout = 3600)
#' temp.dir <- file.path(tempdir(), "HIC_DATA")
#' dir.create(temp.dir)
#' Hic.url <- paste0("https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/",
#'   "7386f953-8da9-47b0-acb2-931cba810544/4DNFIOTPSS3L.hic")
#' HicOutput.pth <- file.path(temp.dir, "Control_HIC.hic")
#' download.file(Hic.url, HicOutput.pth, method = 'curl', extra = '-k')
#' # Import file in R
#' HiC.cmx_lst <- ImportHiC(
#'   file.pth=HicOutput.pth,
#'   res.num=1000,
#'   chrom_1.chr=c("2L", "2R", "2L"),
#'   chrom_2.chr=c("2L", "2R", "2L")
#' )
#' }
ImportHiC <- function(file.pth=NULL, res.num=NULL, chromSize.dtf=NULL, chrom_1.chr=NULL, chrom_2.chr=NULL, verbose.bln=TRUE, cores.num=1){
    # Resolution Format
        options(scipen=999)
        if(inherits(res.num,"character")){
            res.num <- GenomicSystem(res.num)
        }
    # Chromosomes Format
        if(is.null(chrom_2.chr)){
            chrom_2.chr <- chrom_1.chr
        }else if (length(chrom_1.chr) != length(chrom_2.chr)){
            stop("ERROR: chrom_1.chr and chrom_2.chr must have the same length")
        }
        chrom.chr <- c(chrom_1.chr,chrom_2.chr) |> unlist() |> unique()
        if(grepl(pattern="chr", chrom.chr[1], fixed=TRUE)){
            seqlevelsStyleHiC <- "UCSC"
        }else{
            seqlevelsStyleHiC <- "ensembl"
        }
    # Get SeqInfo
        if(GetFileExtension(file.pth)=="hic"){ 
                if("index" %in% colnames(chromSize.dtf)){
                    chromSize.dtf <- strawr::readHicChroms(file.pth) |> dplyr::select(-"index")
                }else{
                    chromSize.dtf <- strawr::readHicChroms(file.pth)
                }
        }else if(GetFileExtension(file.pth) %in% c("cool","mcool", "HDF5", "hdf5", "h5")){ 
            # Define HDF5groups
                chr.group <- ifelse(GetFileExtension(file.pth) %in% c("cool", "HDF5", "hdf5", "h5"),yes="/chroms",no=paste('resolutions',res.num,'chroms',sep='/')) 
            # Get SeqInfo
                chromSize.dtf <- data.frame(rhdf5::h5read(file.pth,name=chr.group))
                rownames(chromSize.dtf) <- chromSize.dtf$name
        }else if(GetFileExtension(file.pth)=="bedpe" & !is.null(chromSize.dtf)){ 
            hic.gnp <- rtracklayer::import(file.pth, format="bedpe")
            megaHic.dtf <- data.frame(
                chrom_1 = as.vector(hic.gnp@first@seqnames),
                i=ceiling(hic.gnp@first@ranges@start/res.num),
                chrom_2 = as.vector(hic.gnp@second@seqnames),
                j=ceiling(hic.gnp@second@ranges@start/res.num),
                counts=hic.gnp@elementMetadata$score
            )
        }else{
            stop("ERROR: file must be .hic, .cool, .mcool, .hdf5, .HDF5 or .bedpe")
        }
        # Standardize seqlevelsStyle of chromSize.dtf according to chrom.chr
            if(grepl("chr", rownames(chromSize.dtf)[1], fixed=TRUE) & seqlevelsStyleHiC=="ensembl"){
                chromSize.dtf$name <- unlist(lapply(strsplit(rownames(chromSize.dtf), "chr"), `[[`, 2))
                rownames(chromSize.dtf) <- unlist(lapply(strsplit(rownames(chromSize.dtf), "chr"), `[[`, 2))
            }else if(!grepl("chr", rownames(chromSize.dtf)[1], fixed=TRUE) & seqlevelsStyleHiC=="UCSC"){
                chromSize.dtf$name <- paste0("chr", rownames(chromSize.dtf))
                rownames(chromSize.dtf) <- paste0("chr", rownames(chromSize.dtf))
            }
        # Standardize chromSize.dtf and chrom.chr
            chromSize.dtf <- dplyr::filter(chromSize.dtf, chromSize.dtf$name %in% chrom.chr)
            chromSize.dtf <- dplyr::mutate(chromSize.dtf, "dimension"=ceiling(chromSize.dtf$length/res.num))
            chrom.chr <- chrom.chr[chrom.chr %in% chromSize.dtf$name]
            chrom_1.chr <- chrom_1.chr[chrom_1.chr %in% chromSize.dtf$name]
            chrom_2.chr <- chrom_2.chr[chrom_2.chr %in% chromSize.dtf$name]
    # Create Genome as GRanges
        binnedGenome.grn <- chromSize.dtf |>
            dplyr::pull("length") |>
            stats::setNames(chromSize.dtf$name) |>
            GenomicRanges::tileGenome(tilewidth=res.num, cut.last.tile.in.chrom=TRUE)
        GenomeInfoDb::seqlengths(binnedGenome.grn) <- chromSize.dtf$length |> stats::setNames(chromSize.dtf$name) 
            chromComb.lst <- paste(chrom_1.chr,chrom_2.chr, sep="_")
        matrixSymmetric.bln <- strsplit(chromComb.lst,"_") |>
            lapply(function(name.chr){
                name.chr[[1]]  == name.chr[[2]]
                }
            ) |>
            unlist()
        matrixType.str <- chromComb.lst
        matrixType.str[which(matrixSymmetric.bln)] <- "cis"
        matrixType.str[which(!matrixSymmetric.bln)] <- "trans"
        matrixKind.str <- chromComb.lst
        matrixKind.str[which(matrixSymmetric.bln)] <- "U"
        matrixKind.str[which(!matrixSymmetric.bln)] <- NA
        attributes.tbl <- dplyr::bind_cols(
            name=chromComb.lst,
            type=matrixType.str,
            kind=matrixKind.str,
            symmetric=matrixSymmetric.bln
        )
    # Dump file
        multicoreParam <- MakeParallelParam(cores.num = cores.num, verbose.bln = verbose.bln)
        hic.lst_cmx <- BiocParallel::bplapply(BPPARAM = multicoreParam, seq_along(chromComb.lst), function(ele.ndx){
            # Chromosomes
                ele.lst <- unlist(strsplit(chromComb.lst[[ele.ndx]],"_"))
                chrom_1.chr <- ele.lst[[1]]
                chrom_2.chr <- ele.lst[[2]]
            # Dimension   
                dims.num <- ele.lst |>
                    lapply(function(chrom){
                        dplyr::filter(chromSize.dtf, chromSize.dtf$name ==chrom) |>
                            dplyr::pull("dimension")
                        }) |>
                    unlist()
            if(GetFileExtension(file.pth)=="hic"){ 
                # Read .hic file
                    hic.dtf <- strawr::straw("NONE", file.pth, chrom_1.chr,chrom_2.chr, "BP", res.num, "observed")
                    hic.dtf$j <- ceiling((hic.dtf$y+1) / res.num)
                    hic.dtf$i <- ceiling((hic.dtf$x+1) / res.num)
            }else if(GetFileExtension(file.pth) %in% c("cool","mcool", "HDF5", "hdf5", "h5")){
                # Define HDF5groups
                    indexes.group <- ifelse(GetFileExtension(file.pth) %in% c("cool", "HDF5", "hdf5", "h5"),yes="/indexes",no=paste('resolutions',res.num,'indexes',sep='/')) 
                    pixels.group <- ifelse(GetFileExtension(file.pth) %in% c("cool", "HDF5", "hdf5", "h5"),yes="/pixels",no=paste('resolutions',res.num,'pixels',sep='/')) 
                # Define start and end of chromosomes
                    ends.ndx <- chromSize.dtf$dimension |> cumsum() |> stats::setNames(chromSize.dtf$name)
                    starts.ndx <- 1+c(0,ends.ndx[-length(ends.ndx)]) |> stats::setNames(chromSize.dtf$name)
                # Read .mcool file
                    bin1.ndx <- as.vector(rhdf5::h5read(file.pth,name=paste(indexes.group,'bin1_offset',sep="/"),index=list(starts.ndx[chrom_1.chr]:ends.ndx[chrom_1.chr])))
                    slice.num <- sum(bin1.ndx[-1] - bin1.ndx[-length(bin1.ndx)])-1
                    chunk.num  <- seq(bin1.ndx[1]+1,bin1.ndx[1]+1+slice.num)
                    hic.dtf <- data.frame(
                        i = as.vector(rhdf5::h5read(file.pth,
                            name=paste(pixels.group,'bin1_id',sep="/"),
                            index=list(chunk.num)))+1,
                        j = as.vector(rhdf5::h5read(file.pth,
                            name=paste(pixels.group,'bin2_id',sep="/"),
                            index=list(chunk.num)))+1,
                        counts = as.vector(rhdf5::h5read(file.pth,
                            name=paste(pixels.group,'count',sep="/"),
                            index=list(chunk.num)))
                    )
                    filter.bin2 <- hic.dtf$j %in% starts.ndx[chrom_2.chr]:ends.ndx[chrom_2.chr]
                    hic.dtf <- hic.dtf[filter.bin2,]
                    hic.dtf <- dplyr::mutate(hic.dtf, i=hic.dtf$i-starts.ndx[chrom_1.chr]+1)
                    hic.dtf <- dplyr::mutate(hic.dtf, j=hic.dtf$j-starts.ndx[chrom_2.chr]+1)         
            }else if(GetFileExtension(file.pth) == "bedpe"){ 
                hic.dtf <-  dplyr::filter(megaHic.dtf, megaHic.dtf$chrom_1==chrom_1.chr & megaHic.dtf$chrom_2==chrom_2.chr)
            }
            # Create Contact matrix
                hic.spm <- Matrix::sparseMatrix(i=hic.dtf$i, j=hic.dtf$j, x=hic.dtf$counts, dims=dims.num )
                row.regions = binnedGenome.grn[which(as.vector(binnedGenome.grn@seqnames) == chrom_1.chr)]
                col.regions = binnedGenome.grn[which(as.vector(binnedGenome.grn@seqnames) == chrom_2.chr)]
                hic.cmx = InteractionSet::ContactMatrix(hic.spm, row.regions, col.regions)
            # Metadata
                hic.cmx@metadata <- dplyr::filter(attributes.tbl, attributes.tbl$name == paste(ele.lst,collapse="_")) |>
                        tibble::add_column(resolution = res.num) |>
                        as.list()
            return(hic.cmx)
        })
    # Add attributes
        hic.lst_cmx <- hic.lst_cmx |>
            stats::setNames(chromComb.lst) |>
            AddAttr(list(
                resolution = res.num,
                chromSize = tibble::as_tibble(chromSize.dtf),
                matricesKind=attributes.tbl,
                mtx="obs"
                ))
        return(hic.lst_cmx)
}