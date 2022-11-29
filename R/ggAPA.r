#' Aggregation plot
#'
#' ggAPA
#' @description Create a ggplot object used for plot aggregation.
#' @param apa.mtx <matrix> : The matrix to plot. (Default NULL)
#' @param title.chr <character> : The title of plot. (Default NULL)
#' @param trimPrct.num <numeric> : A number between 0 and 100 that give the percentage of trimming. (Default 0)
#' @param bounds.chr <character> : Which boundary must be trim, if it's both, trim half of the percentage in inferior and superior see QtlThreshold. (Default "both")
#' @param minBoundary.num <numeric> : Minimal value of Heatmap, force color range. If Null automaticaly find. (Default NULL)
#' @param center.num <numeric> : Center value of Heatmap, force color range. If Null automaticaly find. (Default NULL)
#' @param maxBoundary.num <numeric> : Maximal value of Heatmap, force color range. If Null automaticaly find. (Default NULL)
#' @param colBreaks.num <numeric> : Repartition of colors. If Null automaticaly find. (Default NULL)
#' @param blurPass.num <numeric> : Number of blur pass. (Default 0)
#' @param blurBox.num <numeric> : if null automaticaly compute for 3 Sd. (Default NULL)
#' @param blurSize.num <numeric> : Size of box applied to blurr if null automaticaly compute for 3 Sd. (Default NULL)
#' @param blurSd.num <numeric> : SD of gaussian smooth. (Default 0.5)
#' @param lowerTri.num <numeric> : The value that replace all value in the lower triangle of matrice (Usefull when blur is apply). (Default NULL)
#' @param heatmap.col <character> : Heatmap color list. If null automaticaly compute. (Default NULL)
#' @param na.col <character> : color of NA values. (Default "#F2F2F2")
#' @param colorScale.chr <character> :  shape of color scale on of "linear" or "density" based. (Default "linear")
#' @param bias.num <numeric> : a positive number.  Higher values give more widely spaced colors at the high end. See ?grDevices::colorRamp for more details. (Default 1)
#' @param paletteLength.num <numeric> : The number of color in the palette. (Default 51)
#' @return A ggplot object.
#' @examples
#' library(GenomicED)
#' data(aggreg.mtx)
#'
#'
#' ggAPA(
#'     apa.mtx      = aggreg.mtx,
#'     title.chr    = "APA center on 0",
#'     center.num   = 0,
#'     trimPrct.num = 5,
#'     bounds.chr   = "both",
#'     blurPass.num = 1,
#'     blurSd.num   = 0.5,
#'     heatmap.col  = viridis(6)
#' )

ggAPA = function(
        apa.mtx = NULL,
        title.chr = NULL,
        trimPrct.num=0, # NULL or [0-100]. Percentage of trimming
        bounds.chr="both", # inferior, bot # which boundary must be trim, if it's both, trim half of the percentage in inferior and superior see QtlThreshold
        minBoundary.num=NULL, # Minimal value of Heatmap, force color range. If Null automaticaly find
        center.num=NULL, # Center value of Heatmap, force color range.  If Null automaticaly find
        maxBoundary.num=NULL, # Maximal value of Heatmap, force color range.  If Null automaticaly find
        colBreaks.num=NULL, # Repartition of colors. If Null automaticaly find
        blurPass.num=0, # Number of blur pass
        blurBox.num=NULL, # if null automaticaly compute for 3 Sd
        blurSize.num=NULL, # Size of box applied to blurr if null automaticaly compute for 3 Sd
        blurSd.num=0.5, # SD of gaussian smooth
        lowerTri.num=NULL,
        heatmap.col=NULL, # Heatmap color list. If null automaticaly compute
        na.col="#F2F2F2", # color of NA values
        colorScale.chr="linear",# or "density" ; shape of color scale
        bias.num=1, # bias of color scale see colorRampPalette
        paletteLength.num = 51
    ){
        #############
        # Trimming
        #############
                if(!is.null(colBreaks.num)){
                    minBoundary.num <- min(colBreaks.num)
                    maxBoundary.num <- max(colBreaks.num)
                }
                vec.num <- c(apa.mtx)
                if(is.null(trimPrct.num)){trimPrct.num <- 0}
                if(trimPrct.num!=0 || !is.null(minBoundary.num) || !is.null(maxBoundary.num)){
                    bounds.num_vec <- vec.num |>
                        QtlThreshold(prct.num=trimPrct.num, bounds.chr=bounds.chr) |>
                        stats::setNames(NULL)
                    bounds.num_lst <- list(bounds.num_vec, list(minBoundary.num, maxBoundary.num))
                    bounds.num_lst <- TransposeList(bounds.num_lst)
                    bounds.num_vec <- c(max(unlist(bounds.num_lst[1]),na.rm=TRUE),min(unlist(bounds.num_lst[2]), na.rm=TRUE))
                }else{
                    bounds.num_vec <- NULL
                }
                if(!is.null(bounds.num_vec)){
                    apa.mtx <- TrimOutliers(x.num=apa.mtx,tresholds.num=bounds.num_vec, clip.bln=TRUE)
                    vec.num <- c(apa.mtx)
                }
        #############
        # Smoothing
        #############
            if (blurPass.num){
                for(i in seq_len(blurPass.num)){
                    apa.mtx <- BoxBlur(mat.mtx=apa.mtx, sd.num=blurSd.num, box.num=blurBox.num, boxSize.num=blurSize.num)
                }
                if(!is.null(lowerTri.num)){
                    apa.mtx[lower.tri(apa.mtx, diag=FALSE)] <- lowerTri.num
                }
                vec.num <- c(apa.mtx)
            }
        #############
        # Breaks
        #############
            if(is.null(colBreaks.num)){
                colBreaks.num <- BreakVector(
                    x.num=vec.num,
                    min.num=minBoundary.num,
                    center.num=center.num,
                    max.num=maxBoundary.num,
                    n.num=paletteLength.num,
                    method.chr=colorScale.chr
                )
                minBoundary.num <- min(colBreaks.num)
                maxBoundary.num <- max(colBreaks.num)
            }
        #############
        # Colors
        #############
            if(is.null(heatmap.col)){
                heatmap.col <- dplyr::case_when(
                    !is.null(center.num) && max(colBreaks.num)<=center.num  ~ rev(YlGnBu(paletteLength.num=paletteLength.num,bias=bias.num)),
                    !is.null(center.num) && center.num<=min(colBreaks.num)  ~ YlOrRd(paletteLength.num=paletteLength.num,bias=bias.num ),
                    TRUE                                                    ~ c(rev(YlGnBu(paletteLength.num=floor((paletteLength.num-1)/2),bias=bias.num)), "#FFFFD8", YlOrRd(paletteLength.num=ceiling((paletteLength.num-1)/2), bias=bias.num))
                )
            }
        #############
        # Raster
        #############
            data.dtf <- MeltSpm(apa.mtx)
            plot.ggp <- ggplot2::ggplot(data.dtf, ggplot2::aes(data.dtf$j, data.dtf$i)) +
                ggplot2::geom_raster(ggplot2::aes(fill=data.dtf$x)) +
                ggplot2::scale_fill_gradientn(colours=heatmap.col, values=MinMaxScale(colBreaks.num), na.value=na.col, limits=c(minBoundary.num,maxBoundary.num)) +
                ggplot2::scale_y_reverse(breaks=seq_along(colnames(apa.mtx)), labels=colnames(apa.mtx)) +
                ggplot2::scale_x_continuous(breaks=seq_along(rownames(apa.mtx)), labels=rownames(apa.mtx)) +
                ggplot2::labs(title=title.chr, y=dimnames(apa.mtx)[[2]], x=dimnames(apa.mtx)[[2]]) +
                ggplot2::theme_classic() + ggplot2::theme(
                    axis.line.y=ggplot2::element_blank(), axis.ticks.y=ggplot2::element_blank(),
                    axis.line.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_blank(),
                    legend.title=ggplot2::element_blank()
                )
            return(plot.ggp)
}