#' Configure parallel parameters.
#'
#' MakeParallelParam
#' @description Create BiocParallel parameter according to OS.
#' @param cores.num <numerical> : An integer to specify the number of cores. (Default 1)
#' @param verbose.bln <logical>: A logical value. If TRUE show the progression in console. (Default TRUE)
#' @return  return parrallel parameter according number of cores and OS to use with BiocParallel package.
#' @examples
#' multicoreParam <- MakeParallelParam(2)
#' BiocParallel::bplapply(BPPARAM = multicoreParam, seq_len(3), sqrt)

MakeParallelParam <- function(cores.num = 1, verbose.bln = FALSE){
    if(!is.numeric(cores.num) | cores.num<2 | .Platform$OS.type=="windows"){
        return(BiocParallel::SerialParam(progressbar = verbose.bln))
    }else{
        return(BiocParallel::MulticoreParam(workers = cores.num, progressbar = verbose.bln))
    }   
}