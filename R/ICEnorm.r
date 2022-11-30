#' Compute Iterative Correction on hic chunk.
#' 
#' ICEnorm
#' @description Compute Iterative Correction (Vanilla Count) on hic maps.
#' @param hic.cmx <contactMatrix>: The HiC maps chunk to normalize.
#' @param qtlTh.num <numerical>: the threshold quantile below which the bins will be ignored. (Default 0.15)
#' @param maxIter.num <numerical>: The maximum iteration number.
#' @return A normalized contactMatrix
#' @examples
#' data(HiC_Ctrl.cmx_lst)
#' HiC_Ctrl_Norm.cmx <- ICEnorm(HiC_Ctrl.cmx_lst[['2L_2L']])

ICEnorm <- function(hic.cmx,  qtlTh.num=0.15, maxIter.num=50){
    # Removed Low counts bins
        if(qtlTh.num){
            if(hic.cmx@metadata$symmetric){
                rowBias.num <- Matrix::rowSums(hic.cmx@matrix,na.rm=TRUE) + Matrix::colSums(hic.cmx@matrix,na.rm=TRUE) - Matrix::diag(hic.cmx@matrix)
                colBias.num <- rowBias.num
            }else{
                rowBias.num <- Matrix::rowSums(hic.cmx@matrix,na.rm=TRUE)
                colBias.num <- Matrix::colSums(hic.cmx@matrix,na.rm=TRUE)
            }
            row.ndx <- which(rowBias.num<stats::quantile(rowBias.num,qtlTh.num) & rowBias.num > 0)
            col.ndx <- which(colBias.num<stats::quantile(colBias.num,qtlTh.num) & colBias.num > 0)
            
            meltedHic.dtf <- MeltSpm(hic.cmx@matrix)

            removedHic.dtf <- dplyr::filter(meltedHic.dtf, meltedHic.dtf$i %in% row.ndx | meltedHic.dtf$j %in% col.ndx)
            hic.cmx@metadata$removedCounts <- Matrix::sparseMatrix(i=removedHic.dtf$i,j=removedHic.dtf$j,x=removedHic.dtf$x,dims=dim(hic.cmx@matrix))
        
            hic.dtf <- dplyr::filter(meltedHic.dtf, NotIn(meltedHic.dtf$i,row.ndx) & NotIn(meltedHic.dtf$j,col.ndx))
            hic.cmx@matrix <- Matrix::sparseMatrix(i=hic.dtf$i,j=hic.dtf$j,x=hic.dtf$x,dims=dim(hic.cmx@matrix))
        }
    observed.num <- hic.cmx@matrix@x
    bias.num <- Matrix::rowSums(hic.cmx@matrix,na.rm=TRUE) + Matrix::colSums(hic.cmx@matrix,na.rm=TRUE) - Matrix::diag(hic.cmx@matrix)
    iter.num <- 1
    fit.lm <- stats::lm(c(stats::var(bias.num[which(bias.num!=0)]), 0) ~ c(iter.num, maxIter.num))
    slope.num = stats::coef(fit.lm)[2]
    angle.num = 90-tan(slope.num) * pi /180
    intercept.num = stats::coef(fit.lm)[1]
    max_gain.num = -Inf
    while(iter.num < maxIter.num){
        hic.cmx <- VCnorm(hic.cmx, qtlTh.num=0, sqrt.bln=FALSE)
        bias.num  <- Matrix::rowSums(hic.cmx@matrix,na.rm=TRUE) + Matrix::colSums(hic.cmx@matrix,na.rm=TRUE) - Matrix::diag(hic.cmx@matrix)
        vertical_distances.num = abs(slope.num*iter.num+intercept.num - stats::var(bias.num[which(bias.num !=0)]))
        perpendicular_distance.num = sin(angle.num) * vertical_distances.num
        gain.num = vertical_distances.num - perpendicular_distance.num
        if(gain.num>max_gain.num){
            max_gain.num <- gain.num
            i_max <- iter.num
            hicNorm.cmx <- hic.cmx
        }else if (i_max+floor(maxIter.num*0.1) <= iter.num ){
            break
        }
        iter.num <- iter.num + 1
    }
    hicNorm.cmx@metadata$normalizer <- hicNorm.cmx@matrix@x / observed.num
    hicNorm.cmx@metadata$observed <- observed.num
    return(hicNorm.cmx)
}