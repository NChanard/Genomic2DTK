#' Comoute HiC map normalization.
#' 
#' NormalizeHiC
#' @description Apply a normalization method to a list of contacts matrix.
#' @param hic.cmx_lst <List[contactMatrix]>: The HiC maps list.
#' @param method.chr <character> : The kind of normalization method. One of "ICE", "VC" or "VC_SQRT" (Default "ICE")
#' @param type.chr <character> : Whether normalization must be apply on "cis" chunks, "trans" chunks or "all" chuncks. (Default "cis")
#' @param maxIter.num <numerical>: The maximum iteration number.
#' @param qtlTh.num <numerical>: the threshold quantile below which the bins will be ignored. (Default 0.15)
#' @param cores.num <numerical> : An integer to specify the number of cores. (Default 1)
#' @param verbose.bln <logical>: A logical value. If TRUE show the progression in console. (Default TRUE)
#' @return A matrices list.
#' @examples
NormalizeHiC <- function(hic.cmx_lst, method.chr="ICE", type.chr=c("cis", "trans", "all"), maxIter.num=50, qtlTh.num=0.15, cores.num=1, verbose.bln=TRUE ){
    if (type.chr=="all"){
        megaHic.cmx <- JoinHiC(hic.cmx_lst)
        if(method.chr=="VC"){
            megaHic.cmx <-  VCnorm(megaHic.cmx, qtlTh.num=qtlTh.num, sqrt.bln=FALSE)
        }else if(method.chr=="VC_SQRT"){
            megaHic.cmx <-  VCnorm(megaHic.cmx, qtlTh.num=qtlTh.num, sqrt.bln=TRUE)
        }else if (method.chr=="ICE"){
            megaHic.cmx <- ICEnorm(megaHic.cmx, qtlTh.num=qtlTh.num, maxIter.num=maxIter.num)
        }
        hic.cmx_lst <- CutHiC(megaHic.cmx, verbose.bln=verbose.bln)
    }else if(type.chr=="cis"){
        matricesKind.tbl <- attributes(hic.cmx_lst)$matricesKind
        matricesKind.tbl <- dplyr::filter(matricesKind.tbl, matricesKind.tbl$type=="cis")
        matricesNames.chr <- matricesKind.tbl |> dplyr::pull("name")
        if(cores.num==1){
            start.tim <- Sys.time()
            jobLenght.num <- length(matricesNames.chr)
            hic.cmx_lst[matricesNames.chr] <- lapply(seq_along(matricesNames.chr), function(ele.ndx){
                if(verbose.bln){SuperTK::ShowLoading(start.tim, ele.ndx,jobLenght.num)}
                matrixName.chr <- matricesNames.chr[[ele.ndx]]
                if(method.chr=="VC"){
                    hic.cmx <-  VCnorm(hic.cmx_lst[[matrixName.chr]], qtlTh.num=qtlTh.num, sqrt.bln=FALSE)
                }else if(method.chr=="VC_SQRT"){
                    hic.cmx <-  VCnorm(hic.cmx_lst[[matrixName.chr]], qtlTh.num=qtlTh.num, sqrt.bln=TRUE)
                }else if (method.chr=="ICE"){
                    hic.cmx <- ICEnorm(hic.cmx_lst[[matrixName.chr]], qtlTh.num=qtlTh.num, maxIter.num=maxIter.num)
                }
                return(hic.cmx)
            })
        }else if (cores.num>=2){
            parCl <- parallel::makeCluster(cores.num, type ="FORK")
            hic.cmx_lst[matricesNames.chr] <- parallel::parLapply(parCl,seq_along(matricesNames.chr), function(ele.ndx){
                matrixName.chr <- matricesNames.chr[[ele.ndx]]
                if(method.chr=="VC"){
                    hic.cmx <-  VCnorm(hic.cmx_lst[[matrixName.chr]], qtlTh.num=qtlTh.num, sqrt.bln=FALSE)
                }else if(method.chr=="VC_SQRT"){
                    hic.cmx <-  VCnorm(hic.cmx_lst[[matrixName.chr]], qtlTh.num=qtlTh.num, sqrt.bln=TRUE)
                }else if (method.chr=="ICE"){
                    hic.cmx <- ICEnorm(hic.cmx_lst[[matrixName.chr]], qtlTh.num=qtlTh.num, maxIter.num=maxIter.num)
                }
                return(hic.cmx)
            })
            parallel::stopCluster(parCl)
        }
    }else if(type.chr=="trans"){
        matricesKind.tbl <- attributes(hic.cmx_lst)$matricesKind
        matricesKind.tbl <- dplyr::filter(matricesKind.tbl, matricesKind.tbl$type=="trans")
        matricesNames.chr <- matricesKind.tbl |> dplyr::pull("name")
        chromNames.chr <- matricesNames.chr |> strsplit("_") |> unlist() |> unique()
        chromSize.tbl <- attributes(hic.cmx_lst)$chromSize
        chromSize.tbl <- dplyr::filter(chromSize.tbl, chromSize.tbl$name %in% chromNames.chr)
        trans.cmx_lst <- hic.cmx_lst[matricesNames.chr] |> SuperTK::AddAttr(list(
            resolution   = attributes(hic.cmx_lst)$resolution,
            chromSize    = chromSize.tbl,
            matricesKind = matricesKind.tbl
        ))
        megaHic.cmx <- JoinHiC(trans.cmx_lst)
        if(method.chr=="VC"){
            megaHic.cmx <-  VCnorm(megaHic.cmx, qtlTh.num=qtlTh.num, sqrt.bln=FALSE)
        }else if(method.chr=="VC_SQRT"){
            megaHic.cmx <-  VCnorm(megaHic.cmx, qtlTh.num=qtlTh.num, sqrt.bln=TRUE)
        }else if (method.chr=="ICE"){
            megaHic.cmx <- ICEnorm(megaHic.cmx, qtlTh.num=qtlTh.num, maxIter.num=maxIter.num)
        }
        trans.cmx_lst <- CutHiC(megaHic.cmx, verbose.bln=verbose.bln)
        hic.cmx_lst[matricesNames.chr] <- trans.cmx_lst[matricesNames.chr]
    }
    return(hic.cmx_lst)
}