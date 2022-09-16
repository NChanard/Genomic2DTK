#' FilterInteractions
#'
#' Search in a GInteraction object which interactions correspond ti a target list and return a list of index or filter a matrices list according to target and a selection function.
#' @param matrices.lst <List[matrix]>: The matrices list to filter. If is not NULL, the function will return the filtred matrices list, else return a list of index.
#' @param interarctions.gni <GInteractions>: The GInteraction object on which compute the filter.
#' @param target.lst <List>: a nammed list that describe the target.
#' @param selection.fun <function>: A function that described how the target must be cross. (Defaul intersection of all targets)
#' @return A list of elements index or a filtred matrices list with attributes updates.

FilterInteractions = function(matrices.lst=NULL, interarctions.gni=NULL, target.lst=NULL, selection.fun=function(){Reduce(intersect,interarctions.ndx_lst)}) {
    if(!is.null(matrices.lst) && !is.null(attributes(matrices.lst)$interactions)){interarctions.gni <- attributes(matrices.lst)$interactions}
    interarctions.ndx_lst <- lapply(seq_along(target.lst), function(target.ndx){
        columnName.chr = names(target.lst)[target.ndx]
        if(is.function(target.lst[[target.ndx]])){
            lapply(S4Vectors::mcols(interarctions.gni)[,columnName.chr], target.lst[[target.ndx]]) %>% unlist %>% which
        }else if(is.character(target.lst[[target.ndx]])){
            lapply(S4Vectors::mcols(interarctions.gni)[,columnName.chr],function(columnElement){intersect(as.character(columnElement),target.lst[[columnName.chr]]) %>% length %>% as.logical %>% return(.data)}) %>% unlist %>% which
        }else if(inherits(target.lst[[target.ndx]],"GRanges")){
            GenomicRanges::findOverlaps(InteractionSet::anchors(interarctions.gni)[[columnName.chr]],target.lst[[target.ndx]])@from
        }else if(inherits(target.lst[[target.ndx]],"GInteractions")){
            InteractionSet::findOverlaps(interarctions.gni,target.lst[[target.ndx]])@from
        }
    }) %>% magrittr::set_names(names(target.lst))
    attach(interarctions.ndx_lst)
    on.exit(detach(interarctions.ndx_lst))
    if(length(target.lst)==1){
        interarctions.ndx <- unlist(interarctions.ndx_lst)
    }else if(!is.null(selection.fun)){
        interarctions.ndx <- selection.fun()
    }else{
        return(interarctions.ndx_lst)
    }
    if(!is.null(matrices.lst)){
        matrices.lst[interarctions.ndx] %>%
            DevTK::AddAttr(attribute.lst=list(interactions=interarctions.gni[interarctions.ndx], target=target.lst, selection=selection.fun)) %>%
            DevTK::AddAttr(attribute.lst=attributes(matrices.lst)) %>%
            return(.data)
    }else{
        return(interarctions.ndx)
    }
}