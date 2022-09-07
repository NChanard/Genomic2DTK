#' ggAPA
#'
#' Create a ggplot object used for plot aggragation.
#' @param apa.mtx <matrix> : The matrix to plot. (Default NULL)
#' @param title.chr <character> : The title of plot. (Default NULL)
#' @param trimPrct.num <numeric> : A number between 0 and 100 that give the percentage of trimming. (Default 0)
#' @param bounds.chr <character> : Which boundary must be trim, if it's both, trim half of the percentage in inferior and superior see StatTK::QtlThreshold. (Default "both")
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
#' @param bias.num <numeric> : bias of color scale see grDevices::colorRampPalette (Default 0.8)
#' @param paletteLength.nmb <numeric> : The number of color in the palette. (Default 51)
#' @return A ggplot object.

ggAPA = function(
        apa.mtx = NULL,
        title.chr = NULL,
        trimPrct.num=0, # NULL or [0-100]. Percentage of trimming
        bounds.chr="both", # inferior, bot # which boundary must be trim, if it's both, trim half of the percentage in inferior and superior see StatTK::QtlThreshold
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
        bias.num=0.8, # bias of color scale see colorRampPalette
        paletteLength.nmb = 51
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
                    bounds.num_vec <- vec.num %>%
                        StatTK::QtlThreshold(., prct.num=trimPrct.num, bounds.chr=bounds.chr) %>%
                        magrittr::set_names(NULL) %>%
                        list(., list(minBoundary.num, maxBoundary.num)) %>%
                        purrr::transpose(.) %>% {c(.[1] %>%
                        unlist %>%
                        max(., na.rm=TRUE), .[2] %>%
                        unlist %>% min(., na.rm=TRUE))} %>%
                        suppressWarnings
                }else{
                    bounds.num_vec <- NULL
                }
                if(!is.null(bounds.num_vec)){
                    apa.mtx <- StatTK::TrimOutliers(x.num=apa.mtx,tresholds.num=bounds.num_vec, clip.bln=TRUE)
                    vec.num <- c(apa.mtx)
                }
        #############
        # Smoothing
        #############
            if (blurPass.num){
                for(i in seq_len(blurPass.num)){
                    apa.mtx <- StatTK::BoxBlur(mat.mtx=apa.mtx, sd=blurSd.num, box.num=blurBox.num, boxSize.num=blurSize.num)
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
                colBreaks.num <- StatTK::BreakVector(
                    x.num=vec.num,
                    min.num=minBoundary.num,
                    center.num=center.num,
                    max.num=maxBoundary.num,
                    n.num=paletteLength.nmb,
                    method.chr=colorScale.chr
                )
                minBoundary.num <- min(colBreaks.num)
                maxBoundary.num <- max(colBreaks.num)
            }
        #############
        # Colors
        #############
            if(is.null(heatmap.col)){
                blues.pal <- grDevices::colorRampPalette(colors=rev(RColorBrewer::brewer.pal(9, 'YlGnBu')), space='Lab',interpolate='spline',bias=bias.num)
                reds.pal <- grDevices::colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), space='Lab',interpolate='spline',bias=bias.num)
                heatmap.col <- dplyr::case_when(
                    !is.null(center.num) && max(colBreaks.num)<=center.num  ~ blues.pal(paletteLength.nmb),
                    !is.null(center.num) && center.num<=min(colBreaks.num)  ~ reds.pal(paletteLength.nmb),
                    TRUE                                                    ~ c(blues.pal(floor((paletteLength.nmb-1)/2)),"#FFFFD8",reds.pal(ceiling((paletteLength.nmb-1)/2)))
                )
            }else if(length(heatmap.col)!=paletteLength.nmb){
                heatmap.col <- grDevices::colorRampPalette(colors=heatmap.col, space='Lab', interpolate='spline', bias=bias.num)(paletteLength.nmb)    
            }
        #############
        # Raster
        #############
            ggplot2::ggplot(DataTK::MeltSpm(apa.mtx), ggplot2::aes(j, i)) +
                ggplot2::geom_raster(ggplot2::aes(fill=x)) + 
                ggplot2::scale_fill_gradientn(colours=heatmap.col, values=StatTK::MinMaxScale(colBreaks.num), na.value=na.col, limits=c(minBoundary.num,maxBoundary.num)) +
                ggplot2::scale_y_reverse(breaks=seq_along(colnames(apa.mtx)), labels=colnames(apa.mtx)) +
                ggplot2::scale_x_continuous(breaks=seq_along(rownames(apa.mtx)), labels=rownames(apa.mtx)) +
                ggplot2::labs(title=title.chr, y=dimnames(apa.mtx)[[2]], x=dimnames(apa.mtx)[[2]]) +
                ggplot2::theme_classic() + ggplot2::theme(
                    axis.line.y=ggplot2::element_blank(), axis.ticks.y=ggplot2::element_blank(), 
                    axis.line.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_blank(),
                    legend.title=ggplot2::element_blank()
                ) %>%
            return(.)
}