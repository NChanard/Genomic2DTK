#' Comoute HiC map normalization.
#' 
#' NormalizeHiC
#' @description Apply a normalization method to a list of contacts matrix.
#' @param hic.cmx_lst <List[contactMatrix]>: The HiC maps list.
#' @param method.chr <character> : The kind of normalization method. One of "ICE", "VC" or "VC_SQRT" (Default "ICE")
#' @param type.chr <character> : If NULL normalization is apply on cis chunks then trans chunks. If is "all", normalization is apply on all chuncks at once. (Default NULL)
#' @param maxIter.num <numerical>: The maximum iteration number.
#' @param qtlTh.num <numerical>: the threshold quantile below which the bins will be ignored. (Default 0.15)
#' @param cores.num <numerical> : An integer to specify the number of cores. (Default 1)
#' @param verbose.bln <logical>: A logical value. If TRUE show the progression in console. (Default TRUE)
#' @return A matrices list.
#' @examples
NormalizeHiC <- function(hic.cmx_lst, method.chr="ICE", type.chr=NULL, maxIter.num=50, qtlTh.num=0.15, cores.num=1, verbose.bln=TRUE ){
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
    }else{
        matricesKind.tbl <- attributes(hic.cmx_lst)$matricesKind
        cisMatricesNames.chr <- dplyr::filter(matricesKind.tbl, matricesKind.tbl$type=="cis") |> dplyr::pull("name")
        if(length(cisMatricesNames.chr)){
            if(cores.num==1){
                start.tim <- Sys.time()
                jobLenght.num <- length(cisMatricesNames.chr)
                hic.cmx_lst[cisMatricesNames.chr] <- lapply(seq_along(cisMatricesNames.chr), function(ele.ndx){
                    if(verbose.bln){SuperTK::ShowLoading(start.tim, ele.ndx,jobLenght.num)}
                    matrixName.chr <- cisMatricesNames.chr[[ele.ndx]]
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
                hic.cmx_lst[cisMatricesNames.chr] <- parallel::parLapply(parCl,seq_along(cisMatricesNames.chr), function(ele.ndx){
                    matrixName.chr <- cisMatricesNames.chr[[ele.ndx]]
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
        }
        transMatricesNames.chr <- dplyr::filter(matricesKind.tbl, matricesKind.tbl$type=="trans") |> dplyr::pull("name")
        chromNames.chr <- transMatricesNames.chr |> strsplit("_") |> unlist() |> unique()
        chromSize.tbl <- attributes(hic.cmx_lst)$chromSize
        chromSize.tbl <- dplyr::filter(chromSize.tbl, chromSize.tbl$name %in% chromNames.chr)
        trans.cmx_lst <- hic.cmx_lst[transMatricesNames.chr] |> SuperTK::AddAttr(list(
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
        hic.cmx_lst[transMatricesNames.chr] <- trans.cmx_lst[transMatricesNames.chr]
    }
    return(hic.cmx_lst)
}