#' Aggregation plot
#'
#' ggAPA
#' @description Create a ggplot object used for plot aggregation.
#' @param apa.mtx <matrix> : The matrix to plot. (Default NULL)
#' @param title.chr <character> : The title of plot. (Default NULL)
#' @param trimPrct.num <numeric> : A number between 0 and 100 that give the percentage of trimming. (Default 0)
#' @param bounds.chr <character> : Which boundary must be trim, if it's both, trim half of the percentage in inferior and superior see QtlThreshold. (Default "both")
#' @param colMin.num <numeric> : Minimal value of Heatmap, force color range. If Null automaticaly find. (Default NULL)
#' @param colMid.num <numeric> : Center value of Heatmap, force color range. If Null automaticaly find. (Default NULL)
#' @param colMax.num <numeric> : Maximal value of Heatmap, force color range. If Null automaticaly find. (Default NULL)
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
#' # Index Beaf32 in TADs domains
#' Beaf32_Index.gnr <- IndexFeatures(
#'     gRange.gnr_lst = list(Beaf=Beaf32_Peaks.gnr), 
#'     chromSize.dtf  = data.frame(seqnames = c('2L', '2R'), seqlengths = c(23513712,25286936)),
#'     binSize.num    = 100000
#' )
#' 
#' # Beaf32 <-> Beaf32 Pairing
#' Beaf_Beaf.gni <- SearchPairs(indexAnchor.gnr = Beaf32_Index.gnr)
#' Beaf_Beaf.gni <- Beaf_Beaf.gni[seq_len(2000)] # subset 2000 first for exemple
#' 
#' # Matrices extractions center on Beaf32 <-> Beaf32 point interaction
#' interactions_PF.mtx_lst  <- ExtractSubmatrix(
#'     feature.gn         = Beaf_Beaf.gni,
#'     hic.cmx_lst        = HiC_Ctrl.cmx_lst,
#'     referencePoint.chr = "pf"
#' )
#' 
#' # Aggregate matrices in one matrix
#' aggreg.mtx <- Aggregation(interactions_PF.mtx_lst)
#' 
#' # Visualization
#' ggAPA(
#'     apa.mtx = aggreg.mtx
#' )
#' 
#' # Add Title
#' ggAPA(
#'     apa.mtx = aggreg.mtx,
#'     title.chr = "APA"
#' )
#' 
#' # Trim values
#' ggAPA(
#'     apa.mtx      = aggreg.mtx,
#'     title.chr    = "APA 30% trimmed on upper side of distribution",
#'     trimPrct.num = 30,
#'     bounds.chr   = "upper"
#' )
#' ggAPA(
#'     apa.mtx      = aggreg.mtx,
#'     title.chr    = "APA 30% trimmed on lower side of distribution",
#'     trimPrct.num = 30,
#'     bounds.chr   = "lower"
#' )
#' ggAPA(
#'     apa.mtx      = aggreg.mtx,
#'     title.chr    = "APA 15% trimmed on each side of distribution",
#'     trimPrct.num = 30,
#'     bounds.chr   = "both"
#' )
#' 
#' # Change Minimal, Central and Maximal Colors scale value
#' ggAPA(
#'     apa.mtx    = aggreg.mtx,
#'     title.chr  = "APA [min 200, center 300, max 600]",
#'     colMin.num = 200,
#'     colMid.num = 300,
#'     colMax.num = 600
#' )
#' 
#' # Change Color
#' ggAPA(
#'     apa.mtx     = aggreg.mtx,
#'     title.chr   = "APA",
#'     heatmap.col = viridis(6),
#'     na.col      = "black"
#' )
#' ggAPA(
#'     apa.mtx     = aggreg.mtx,
#'     title.chr   = "APA",
#'     heatmap.col = c("black", "white"),
#' )
#' 
#' # Change Color distribution
#' ggAPA(
#'     apa.mtx       = aggreg.mtx,
#'     title.chr     = "APA [100,150,200,250,300,350,600]",
#'     colBreaks.num = c(100,150,200,250,300,350,600) # Choosen Breaks
#' )
#' ggAPA(
#'     apa.mtx    = aggreg.mtx,
#'     title.chr  = "APA",
#'     colorScale = "density" # color distribution based on density
#' )
#' ggAPA(
#'     apa.mtx   = aggreg.mtx,
#'     title.chr = "APA",
#'     bias.num  = 2 #(>1 wait on extremums)
#' )
#' ggAPA(
#'     apa.mtx   = aggreg.mtx,
#'     title.chr = "APA",
#'     bias.num  = 0.5 # (<1 wait on center)
#' )
#' 
#' # Apply a Blurr
#' ggAPA(
#'     apa.mtx      = aggreg.mtx,
#'     title.chr    = "APA",
#'     blurPass.num = 1,
#'     blurSd.num   = 0.5
#' )
#' 
#' # ggplot2 object modifications  
#' # Since the function returns a ggplot object, it is possible
#' # to modify it following the ggplot2 grammar.
#' ggAPA(
#'     apa.mtx     = aggreg.mtx,
#'     title.chr   = "APA",
#' ) + 
#' ggplot2::labs(
#'     title    = "New title",
#'     subtitle = "and subtitle"
#' )

ggAPA = function(
        apa.mtx = NULL,
        title.chr = NULL,
        trimPrct.num=0, # NULL or [0-100]. Percentage of trimming
        bounds.chr="both", # inferior, bot # which boundary must be trim, if it's both, trim half of the percentage in inferior and superior see QtlThreshold
        colMin.num=NULL, # Minimal value of Heatmap, force color range. If Null automaticaly find
        colMid.num=NULL, # Center value of Heatmap, force color range.  If Null automaticaly find
        colMax.num=NULL, # Maximal value of Heatmap, force color range.  If Null automaticaly find
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
                    colMin.num <- min(colBreaks.num)
                    colMax.num <- max(colBreaks.num)
                }
                vec.num <- c(apa.mtx)
                if(is.null(trimPrct.num)){trimPrct.num <- 0}
                if(trimPrct.num!=0 || !is.null(colMin.num) || !is.null(colMax.num)){
                    bounds.num_vec <- vec.num |>
                        QtlThreshold(prct.num=trimPrct.num, bounds.chr=bounds.chr) |>
                        stats::setNames(NULL)
                    bounds.num_lst <- list(bounds.num_vec, list(colMin.num, colMax.num))
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
                    min.num=colMin.num,
                    center.num=colMid.num,
                    max.num=colMax.num,
                    n.num=paletteLength.num,
                    method.chr=colorScale.chr
                )
                colMin.num <- min(colBreaks.num)
                colMax.num <- max(colBreaks.num)
            }
        #############
        # Colors
        #############
            if(is.null(heatmap.col)){
                heatmap.col <- dplyr::case_when(
                    !is.null(colMid.num) && max(colBreaks.num)<=colMid.num  ~ rev(YlGnBu(paletteLength.num=paletteLength.num,bias=bias.num)),
                    !is.null(colMid.num) && colMid.num<=min(colBreaks.num)  ~ YlOrRd(paletteLength.num=paletteLength.num,bias=bias.num ),
                    TRUE                                                    ~ c(rev(YlGnBu(paletteLength.num=floor((paletteLength.num-1)/2),bias=bias.num)), "#FFFFD8", YlOrRd(paletteLength.num=ceiling((paletteLength.num-1)/2), bias=bias.num))
                )
            }
        #############
        # Raster
        #############
            data.dtf <- MeltSpm(apa.mtx)
            plot.ggp <- ggplot2::ggplot(data.dtf, ggplot2::aes(data.dtf$j, data.dtf$i)) +
                ggplot2::geom_raster(ggplot2::aes(fill=data.dtf$x)) +
                ggplot2::scale_fill_gradientn(colours=heatmap.col, values=MinMaxScale(colBreaks.num), na.value=na.col, limits=c(colMin.num,colMax.num)) +
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