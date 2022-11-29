#' Blur a matrix.
#'
#' BoxBlur
#' @description Blur a matrix with a one dimensional kernel.
#' @param mat.mtx <matrix>: numerical matrix.
#' @param box.num <numeric>: the numerical vector for kernel. If NULL apply a GaussBox (see 'GaussBox' function) (Default NULL)
#' @param boxSize.num <numeric>: if box.num is NULL, size of kernel for 'GaussBox' function. (Default NULL)
#' @param sd.num <numeric>: if box.num is NULL,  standard deviation parameter for 'GaussBox' function. (Default NULL)
#' @return Blurred matrix.
#' @examples
#' set.seed(981643)
#' mat.mtx <- rnorm(10000,50,10)**3 |> matrix(100,100)
#' heatmap(mat.mtx,Rowv=NA,Colv=NA)
#' heatmap(BoxBlur(mat.mtx),Rowv=NA,Colv=NA)

BoxBlur <- function(mat.mtx, box.num=NULL, boxSize.num=NULL, sd.num=1){
    if(is.null(box.num)){box.num <- GaussBox(sd.num=sd.num,scale.chr="1",boxSize.num=boxSize.num)}
    pad.num <- (length(box.num)-1)/2
    mat.mtx <- PadMtx(mat.mtx=mat.mtx, padSize.num=pad.num, value.num=NULL, side.chr=c('top','bot','right','left'))
    matVsmth.mtx2 <- sapply(((1+pad.num):(dim(mat.mtx)[2]-pad.num)), function(j){(t(mat.mtx[,(j-pad.num):(j+pad.num)]) * box.num) |>
        apply(2,Plus)})
        matHsmth.mtx2 <- t(sapply(((1+pad.num):(dim(matVsmth.mtx2)[1]-pad.num)), function(i){
            (matVsmth.mtx2[(i-pad.num):(i+pad.num),] * box.num) |>
            apply(2,Plus) |>
            t()
            }))
    which.ndx <- which(matHsmth.mtx2==0)
    if(length(which.ndx )){
        matHsmth.mtx2 = Rise0(matHsmth.mtx2,which.ndx=which.ndx)
    }
    return(matHsmth.mtx2)
}