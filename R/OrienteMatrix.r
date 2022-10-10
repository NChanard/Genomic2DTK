#' OrienteMatrix
#'
#' Apply a 180° rotation follow with a transposation on a matrix or on matricies in a list according to the interactions attributes of the list.
#' @param matrice.mtx <matrix or List[matrix]>: matrix or matricies list to oriente 
#' @return Oriented matrix or matricies list
#' @examples
#' library(GenomicED)
#' data("submatrixPF_Ctrl.mtx_lst")
#' submatrixPF_Ctrl_oriented.mtx_lst <- OrienteMatrix(submatrixPF_Ctrl.mtx_lst)
OrienteMatrix <- function(matrice.mtx){
    if(is.list(matrice.mtx) && !is.null(attributes(matrice.mtx)$interactions)){
        orientedMatrice.mtx <- matrice.mtx
        orientedMatrice.mtx[which(!attributes(matrice.mtx)$interactions$orientation)] %<>% lapply(OrienteMatrix) 
        orientedMatrice.mtx %<>% SuperTK::AddAttr(attributes(matrice.mtx))
        attributes(orientedMatrice.mtx)$interactions$orientation <- TRUE
        attributes(orientedMatrice.mtx)$interactions$submatrix.name <- attributes(orientedMatrice.mtx)$interactions$name
        names(orientedMatrice.mtx) <- attributes(orientedMatrice.mtx)$interactions$name
        return(orientedMatrice.mtx)
    }else{
        return(t(apply(as.data.frame(apply(matrice.mtx,1,rev)),1,rev)))
    }
}