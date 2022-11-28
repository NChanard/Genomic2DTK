#' Create parallel parameter
#'
#' MakeParallelParam
#' @keywords internal
#' @description Create BiocParallel parameter according to OS.
#' @param cores.num <numerical> : An integer to specify the number of cores. (Default 1)
#' @return  
MakeParallelParam <- function(cores.num = 1, verbose.bln = FALSE){
    if(!is.numeric(cores.num) | cores.num<2 | .Platform$OS.type=="windows"){
        return(BiocParallel::SerialParam(progressbar = verbose.bln))
    }else{
        return(BiocParallel::MulticoreParam(workers = cores.num, progressbar = verbose.bln))
    }   
}
