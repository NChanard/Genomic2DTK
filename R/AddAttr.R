#' Add list as attributes.
#'
#' AddAttr
#' @description Add list as attributes to any object with or without overwrite.
#' @param var.any <any>: an object to which attributes are to be added.
#' @param attribute.lst <list>: a named list of new attributes.
#' @param overwrite.bln <logical>: whether an overwrite is required on attributes with the same name.(Default FALSE)
#' @return the object with new attributes.
#' @examples
#' x <- 1:10
#' x <- AddAttr(x, list(dim=c(2,5)))
#' x
#' x <- AddAttr(x, list(dim=c(5,2)))
#' x
#' x <- AddAttr(x,attribute.lst,overwrite.bln=TRUE)
#' x

AddAttr <- function(var.any=NULL, attribute.lst=NULL, overwrite.bln=FALSE){
    intersectAttr <- intersect(names(attributes(var.any)),  names(attribute.lst))
    if(overwrite.bln & length(intersectAttr)){
        attribute.lst <- c(attributes(var.any)[which(names(attributes(var.any)) != intersectAttr)] , attribute.lst)
    }else if (length(intersectAttr)){
        attribute.lst <- c(attributes(var.any) , attribute.lst[which(names(attribute.lst) != intersectAttr)])
    }else{
        attribute.lst <- c(attributes(var.any) , attribute.lst)
    }
    attributes(var.any) <- attribute.lst
    return(var.any)
}