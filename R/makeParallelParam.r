#' Create parallel parameter
#'
#' makeParallelParam
#' @keywords internal
#' @description Create BiocParallel parameter according to OS.
#' @param cores.num <numerical> : An integer to specify the number of cores. (Default 1)
#' @return  
makeParallelParam <- function(cores.num = 1, verbose.bln = FALSE){
    if(!is.numeric(cores.num) | cores.num<2){
        return(BiocParallel::SerialParam(progressbar = verbose.bln))
    }else{
        if(.Platform$OS.type=="windows"){
            return(BiocParallel::SnowParam(workers = cores.num, progressbar = verbose.bln))
        }else{
            return(BiocParallel::MulticoreParam(workers = cores.num, progressbar = verbose.bln))
        }    
    }
}
