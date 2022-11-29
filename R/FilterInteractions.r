#' Submatrix or Interactions filtering.
#' 
#' FilterInteractions
#' @description Search in a GInteraction object which interactions correspond ti a target list and return a list of index or filter a matrices list according to target and a selection function.
#' @param matrices.lst <List[matrix]>: The matrices list to filter. If is not NULL, the function will return the filtred matrices list, else return a list of index.
#' @param interarctions.gni <GInteractions>: The GInteraction object on which compute the filter.
#' @param target.lst <List>: a nammed list that describe the target.
#' @param selection.fun <function>: A function that described how the target must be cross. (Defaul intersection of all targets)
#' @return A list of elements index or a filtred matrices list with attributes updates.
#' @examples
#' library(GenomicED)
#' data("submatrixPF_Ctrl.mtx_lst")
#'
#' \dontrun{
#' target.lst <- list(
#'   anchor.Beaf.name = c("Beaf32_8","Beaf32_15"),
#'   bait.Tss.name    = c("FBgn0031214","FBgn0005278"),
#'   name             = c("2L:74_2L:77"),
#'   distance         = function(columnElement){
#'     return(14000==columnElement || columnElement == 3000)
#'     }
#' )
#'
#'
#' # Extraction on InteractionSet
#' FilterInteractions(
#'   interarctions.gni = attributes(submatrixPF_Ctrl.mtx_lst)$interactions,
#'   target.lst        = target.lst,
#'   selection.fun     = NULL
#' ) |> str(max.level=1)
#' 
#' selection.fun = function(){
#'   Reduce(intersect, list(anchor.Beaf.name, bait.Tss.name ,distance) ) |>
#'   setdiff(name)
#' }
#'
#' 
#' # Extraction on matrices list and with selection
#' FilterInteractions(
#'   matrices.lst = submatrixPF_Ctrl.mtx_lst,
#'   target.lst        = target.lst,
#'   selection.fun     = selection.fun
#' ) |> str(max.level=1)
#'
#' 
#' # Extraction on InteractionSet with InteractionsSet
#' target.lst <- list(interactions = attributes(submatrixPF_Ctrl.mtx_lst)$interactions[1:2])
#' FilterInteractions(
#'   interarctions.gni = attributes(submatrixPF_Ctrl.mtx_lst)$interactions,
#'   target.lst        = target.lst,
#'   selection.fun     = NULL
#' )
#'
#' 
#' # Extraction on InteractionSet list and with GRanges
#' gRangesTarget.gnr <- attributes(submatrixPF_Ctrl.mtx_lst)$interactions |> 
#' InteractionSet::anchors()
#' target.lst <- list(
#'     first = 
#'       gRangesTarget.gnr[["first"]][1:2]
#' )
#' FilterInteractions(
#'   interarctions.gni = attributes(submatrixPF_Ctrl.mtx_lst)$interactions,
#'   target.lst        = target.lst,
#'   selection.fun     = NULL
#' )
#' }

FilterInteractions = function(matrices.lst=NULL, interarctions.gni=NULL, target.lst=NULL, selection.fun=function(){Reduce(intersect,interarctions.ndx_lst)}) {
    if(!is.null(matrices.lst) && !is.null(attributes(matrices.lst)$interactions)){interarctions.gni <- attributes(matrices.lst)$interactions}
    interarctions.ndx_lst <- lapply(seq_along(target.lst), function(target.ndx){
        columnName.chr = names(target.lst)[target.ndx]
        if(is.function(target.lst[[target.ndx]])){
            lapply(S4Vectors::mcols(interarctions.gni)[,columnName.chr], target.lst[[target.ndx]]) |> unlist() |> which()
        }else if(is.character(target.lst[[target.ndx]])){
            lapply(S4Vectors::mcols(interarctions.gni)[,columnName.chr],function(columnElement){
                res.lgk <- intersect(as.character(columnElement),target.lst[[columnName.chr]]) |> length() |> as.logical()
                return(res.lgk)
                }) |> unlist() |> which()
        }else if(inherits(target.lst[[target.ndx]],"GRanges")){
            GenomicRanges::findOverlaps(InteractionSet::anchors(interarctions.gni)[[columnName.chr]],target.lst[[target.ndx]])@from
        }else if(inherits(target.lst[[target.ndx]],"GInteractions")){
            InteractionSet::findOverlaps(interarctions.gni,target.lst[[target.ndx]])@from
        }
    }) |> stats::setNames(names(target.lst))
    if(length(target.lst)==1){
        interarctions.ndx <- unlist(interarctions.ndx_lst)
    }else if(!is.null(selection.fun)){
        for(target.ndx in seq_along(interarctions.ndx_lst)){
           assign(names(interarctions.ndx_lst)[target.ndx], interarctions.ndx_lst[[target.ndx]], envir = parent.frame())
        }
        interarctions.ndx <- selection.fun()
    }else{
        return(interarctions.ndx_lst)
    }
    if(!is.null(matrices.lst)){
        matrices.filt.lst <- matrices.lst[interarctions.ndx]
        # recover attributes DD221107 # TODO
        attributes(matrices.filt.lst)$interactions = attributes(matrices.lst)$interactions[interarctions.ndx]
        attributes(matrices.filt.lst)$resolution = attributes(matrices.lst)$resolution
        attributes(matrices.filt.lst)$referencePoint = attributes(matrices.lst)$referencePoint
        attributes(matrices.filt.lst)$matriceDim = attributes(matrices.lst)$matriceDim
        attributes(matrices.filt.lst)$target = target.lst
        attributes(matrices.filt.lst)$selection = selection.fun
        #  |>
        #     AddAttr(attribute.lst=list(interactions=interarctions.gni[interarctions.ndx], target=target.lst, selection=selection.fun)) |>
        #     AddAttr(attribute.lst=attributes(matrices.lst))
            return(matrices.filt.lst)
    }else{
        return(interarctions.ndx)
    }
}