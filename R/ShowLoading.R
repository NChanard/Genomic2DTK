#' Show the loading of a looping operation.
#'
#' ShowLoading
#' @description Show the loading of a looping operation.
#' @param start.tim <POSIXct POSIXt>: A time obtain with 'Sys.time()'.
#' @param operation.ndx <numerical>: The index number of the operation.
#' @param operation.num <numerical>: The total number of operation.
#' @examples
#' start.tim <- Sys.time()
#' for(i in seq_len(10000)){
#'     ShowLoading(start.tim, i , 10000)
#' }

ShowLoading = function(start.tim=NULL, operation.ndx=NULL, operation.num=NULL){
    if (operation.ndx==1){
        cat(format(paste0(" ",format(100*operation.ndx/operation.num,digits=2,nsmall=2)," %"),width=9,justify = "left"),"\r")
    }else{
        restingtime.num <- as.numeric( (difftime(Sys.time(), start.tim, units="secs")   / ( operation.ndx/operation.num) ) - difftime(Sys.time(), start.tim, units="secs"))
        if (restingtime.num>3600) {
            restingtime.num <- restingtime.num/3600
            units.chr <- "hours"
        }else if (restingtime.num>60) {
            restingtime.num <- restingtime.num/60
            units.chr <- "mins"
        }else{
            units.chr <- "secs"
        }
        cat(
            format(paste0(" ",format(100* operation.ndx/operation.num,digits=2,nsmall=2)," %"),width=9,justify = "left"),
            paste0(format(paste0("Resting time: ",format(round(restingtime.num,2),digits=2,nsmall=2)),width=19,justify = "left")," ",units.chr),
            "\r")
    }
}
