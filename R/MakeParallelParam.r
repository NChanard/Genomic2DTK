#' Create parallel parameter
#'
#' MakeParallelParam
#' @description Create BiocParallel parameter according to OS.
#' @param cores.num <numerical> : An integer to specify the number of cores. (Default 1)
#' @return  return parrallel parameter according number of cores and OS to use with BiocParallel package.
#' @exemple

MakeParallelParam <- function(cores.num = 1, verbose.bln = FALSE){
    if(!is.numeric(cores.num) | cores.num<2 | .Platform$OS.type=="windows"){
        return(BiocParallel::SerialParam(progressbar = verbose.bln))
    }else{
        return(BiocParallel::MulticoreParam(workers = cores.num, progressbar = verbose.bln))
    }   
}