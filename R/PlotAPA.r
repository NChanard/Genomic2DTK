#' PlotAPA
#'
#' Draw aggregation plot from aggregation matrices.
#' @param apa.mtx <matrix>: The aggregated matrix.
#' @param trimPrct.num <numeric>: A number between 0 and 100 thaht give the percentage of triming in matrices.
#' @param minBoundary.num <matrix>: The minimal value in color scale. If Null automaticaly find.
#' @param center.num <matrix>: The middle value in color scale. If Null automaticaly find.
#' @param maxBoundary.num <matrix>: The mximal value in color scale. If Null automaticaly find.
#' @param minConditionBoundary.num <matrix>: Avalaible for plotting differantial aggregation. The minimal value in color scale in the classsical aggregation plot. If Null automaticaly find.
#' @param maxConditionBoundary.num <matrix>: Avalaible for plotting differantial aggregation. The maxiaml value in color scale in the classsical aggregation plot. If Null automaticaly find.
#' @return None
#' @examples
#' library(GenomicED)
#' data(aggreg.mtx)
#' PlotAPA(
#'     apa.mtx                  = aggreg.mtx,
#'     trimPrct.num             = 20,
#'     minBoundary.num          = -2,
#'     center.num               = 0,
#'     maxBoundary.num          = 2,
#'     minConditionBoundary.num = 0,
#'     maxConditionBoundary.num = 2
#' )
PlotAPA = function(apa.mtx = NULL, trimPrct.num=0, minBoundary.num=NULL, center.num=NULL, maxBoundary.num=NULL, minConditionBoundary.num=NULL, maxConditionBoundary.num=NULL){
    .ggDensity <- function(data.lst=NULL, colour.col=NULL, mean.bln=TRUE, title.chr=NULL){
        data.lst_tbl <- lapply(seq_along(data.lst),function(element.ndx){
            tibble::tibble(
                value = data.lst[[element.ndx]],
                class = factor(names(data.lst)[[element.ndx]])
            ) %>%
            return(.data)
            }) 
        data.tbl <- dplyr::bind_rows(data.lst_tbl)
        if(is.null(colour.col)){
            colour.col <- SuperTK::Hue(length(data.lst)) %>% magrittr::set_names(names(data.lst))
        }
        plot.gp <- ggplot2::ggplot(data.tbl, ggplot2::aes(x=.data$value, fill=.data$class, colour=.data$class)) + 
            ggplot2::geom_density(alpha=0.1) +
            ggplot2::scale_color_manual(values = colour.col)+
            ggplot2::scale_fill_manual(values = colour.col) +
            ggplot2::labs(title=title.chr) 
        if (mean.bln){
            mu.tbl <- data.tbl %>% dplyr::group_by(.data$class) %>% dplyr::summarise(grp.mean = mean(.data$value))
            plot.gp <- plot.gp + 
                ggplot2::geom_vline(data = mu.tbl, ggplot2::aes(xintercept = .data$grp.mean, colour = .data$class), linetype = "dashed")
        }
        return(plot.gp)
    }
    # Differential or not?
        differential.bln <- !is.null(attributes(apa.mtx)$matrices)
        if(differential.bln){
            heatmap.col = NULL
        }else{
            heatmap.col = SuperTK::viridis(255)
        }
    # Plot
        # Auto Scale
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=apa.mtx, 
                heatmap.col=heatmap.col,
                title.chr=ifelse(differential.bln,
                    yes="Agregation of differential matrices",
                    no="Agregation")
            ) + ggplot2::labs(subtitle="scale (auto), center()")
            print(plot.gp)
        # Auto Scale + Center
            if(!is.null(center.num)){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=apa.mtx, 
                    heatmap.col=heatmap.col,
                    center.num=center.num,
                    title.chr=ifelse(differential.bln,
                        yes="Agregation of differential matrices",
                        no="Agregation")
                ) + ggplot2::labs(subtitle=paste0("scale (auto), center(",center.num,")"))
                print(plot.gp)
            }
        # Trim Scale + Center
            if(!is.null(trimPrct.num) && 0<trimPrct.num){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=apa.mtx, 
                    heatmap.col=heatmap.col,
                    trimPrct.num=trimPrct.num,
                    center.num=center.num,
                    title.chr=ifelse(differential.bln,
                        yes="Agregation of differential matrices",
                        no="Agregation")
                ) + ggplot2::labs(subtitle=paste0("scale (rm ",trimPrct.num,"%), center(",center.num,")"))
                print(plot.gp)
            }
        # MinMax Scale + Center
            if(!is.null(minBoundary.num) || !is.null(maxBoundary.num)){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=apa.mtx, 
                    heatmap.col=heatmap.col,
                    minBoundary.num=minBoundary.num,
                    center.num=center.num,
                    maxBoundary.num=maxBoundary.num,
                    title.chr=ifelse(differential.bln,
                        yes="Agregation of differential matrices",
                        no="Agregation")
                ) + ggplot2::labs(subtitle=paste0("scale (",minBoundary.num,";",maxBoundary.num,"), center(",center.num,")"))
                print(plot.gp)
            }
    if (differential.bln){
        # Pval + Auto Scale
            if(!is.null(attributes(apa.mtx)$matrices$pVal) && sum(!is.na(attributes(apa.mtx)$matrices$pVal))>=3){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$pVal,
                    heatmap.col=SuperTK::YlOrRd(9),
                    title.chr = "-log10(p.values)"
                ) + ggplot2::labs(subtitle="scale (auto), center()")
            }else{
                plot.gp <- ggplot2::ggplot() +
                    ggplot2::theme_void() +
                    ggplot2::annotate("text", x = 1, y = 1,
                        label = "Not enough pval computed to plot a pval matrix (<3) or nothing significant")
            }
            print(plot.gp)
        # FiltPval + Auto Scale + Center
            if(!is.null(attributes(apa.mtx)$matrices$aggDiffPvalFilt) && sum(!is.na(attributes(apa.mtx)$matrices$pVal))>=3){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$aggDiffPvalFilt,
                    heatmap.col=heatmap.col,
                    center.num=center.num,
                    title.chr = paste0("Agregation of differential matrices")
                ) + ggplot2::labs(subtitle=paste0("filtred by p.values, scale (auto), center(",center.num,")"))
            }else{
                plot.gp <- ggplot2::ggplot() +
                    ggplot2::theme_void() +
                    ggplot2::annotate("text", x = 1, y = 1,
                    label = "Not enough pval computed to plot a pval matrix (<3) or nothing significant")
            }
            print(plot.gp)
        # FiltPval + Trim Scale + Center
            if(!is.null(attributes(apa.mtx)$matrices$aggDiffPvalFilt) && sum(!is.na(attributes(apa.mtx)$matrices$pVal))>=3){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$aggDiffPvalFilt,
                    heatmap.col=heatmap.col,
                    trimPrct.num=trimPrct.num,
                    center.num=center.num,
                    title.chr = paste0("Agregation of differential matrices")
                ) + ggplot2::labs(subtitle=paste0("filtred by p.values, scale (rm ",trimPrct.num,"%), center(",center.num,")"))
            }else{
                plot.gp <- ggplot2::ggplot() +
                    ggplot2::theme_void() +
                    ggplot2::annotate("text", x = 1, y = 1,
                    label = "Not enough pval computed to plot a pval matrix (<3) or nothing significant")
            }
            print(plot.gp)
        # Delta + Auto Scale + Center
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=attributes(apa.mtx)$matrices$aggDelta, 
                heatmap.col=heatmap.col,
                center.num=center.num,
                title.chr="Differential of agregated matrices"
            ) + ggplot2::labs(subtitle=paste0("scale (auto), center(",center.num,")"))
            print(plot.gp)
        # Delta + Trim Scale + Center
            if(!is.null(trimPrct.num) && 0<trimPrct.num){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$aggDelta, 
                    heatmap.col=heatmap.col,
                    trimPrct.num=trimPrct.num,
                    center.num=center.num,
                    title.chr="Differential of agregated matrices"
                ) + ggplot2::labs(subtitle=paste0("scale (rm ",trimPrct.num,"%), center(",center.num,")"))
                print(plot.gp)
            }
        # Delta + Auto Scale + Center
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=attributes(apa.mtx)$matrices$aggCorrectedDelta, 
                heatmap.col=heatmap.col,
                center.num=center.num,
                title.chr="Differential of corrected agregated matrices"
            ) + ggplot2::labs(subtitle=paste0("scale (auto), center(",center.num,")"))
            print(plot.gp)
        # Delta + Trim Scale + Center
            if(!is.null(trimPrct.num) && 0<trimPrct.num){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$aggCorrectedDelta, 
                    heatmap.col=heatmap.col,
                    trimPrct.num=trimPrct.num,
                    center.num=center.num,
                    title.chr="Differential of corrected agregated matrices"
                ) + ggplot2::labs(subtitle=paste0("scale (rm ",trimPrct.num,"%), center(",center.num,")"))
                print(plot.gp)
            }
        # Control + Auto Scale
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=attributes(apa.mtx)$matrices$aggCtrl, 
                heatmap.col=SuperTK::viridis(51),
                title.chr="Agregation control"
            ) + ggplot2::labs(subtitle="scale (auto)")
            print(plot.gp)
        # Control + Trim Scale
            if(!is.null(trimPrct.num) && 0<trimPrct.num){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$aggCtrl, 
                    heatmap.col=SuperTK::viridis(51),
                    trimPrct.num=trimPrct.num,
                    title.chr="Agregation control"
                ) + ggplot2::labs(subtitle=paste0("scale (rm ",trimPrct.num,"%)"))
                print(plot.gp)
            }
        # Control + MinMax Scale
            if(!is.null(minConditionBoundary.num) || !is.null(maxConditionBoundary.num)){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$aggCtrl, 
                    heatmap.col=SuperTK::viridis(51),
                    minBoundary.num=minConditionBoundary.num,
                    maxBoundary.num=maxConditionBoundary.num,
                    title.chr="Agregation control"
                ) + ggplot2::labs(subtitle=paste0("scale (",minConditionBoundary.num,";",maxConditionBoundary.num,")"))
                print(plot.gp)
            }
        # Condition + Auto Scale
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=attributes(apa.mtx)$matrices$agg, 
                heatmap.col=SuperTK::viridis(51),
                title.chr="Agregation"
            ) + ggplot2::labs(subtitle="scale (auto)")
            print(plot.gp)
        # Condition + Trim Scale
            if(!is.null(trimPrct.num) && 0<trimPrct.num){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$agg, 
                    heatmap.col=SuperTK::viridis(51),
                    trimPrct.num=trimPrct.num,
                    title.chr="Agregation"
                ) + ggplot2::labs(subtitle=paste0("scale (rm ",trimPrct.num,"%)"))
                print(plot.gp)
            }
        # Condition + MinMax Scale
            if(!is.null(minConditionBoundary.num) || !is.null(maxConditionBoundary.num)){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$agg, 
                    heatmap.col=SuperTK::viridis(51),
                    minBoundary.num=minConditionBoundary.num,
                    maxBoundary.num=maxConditionBoundary.num,
                    title.chr="Agregation"
                ) + ggplot2::labs(subtitle=paste0("scale (",minConditionBoundary.num,";",maxConditionBoundary.num,")"))
                print(plot.gp)
            }
        # Corrected condition + Auto Scale
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=attributes(apa.mtx)$matrices$aggCorrected, 
                heatmap.col=SuperTK::viridis(51),
                title.chr="Agregation corrected"
            ) + ggplot2::labs(subtitle="scale (auto)")
            print(plot.gp)
        # Corrected condition + Trim Scale
            if(!is.null(trimPrct.num) && 0<trimPrct.num){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$aggCorrected, 
                    heatmap.col=SuperTK::viridis(51),
                    trimPrct.num=trimPrct.num,
                    title.chr="Agregation corrected"
                ) + ggplot2::labs(subtitle=paste0("scale (rm ",trimPrct.num,"%)"))
                print(plot.gp)
            }
        # Corrected condition + MinMax Scale
            if(!is.null(minConditionBoundary.num) || !is.null(maxConditionBoundary.num)){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$aggCorrected, 
                    heatmap.col=SuperTK::viridis(51),
                    minBoundary.num=minConditionBoundary.num,
                    maxBoundary.num=maxConditionBoundary.num,
                    title.chr="Agregation corrected"
                ) + ggplot2::labs(subtitle=paste0("scale (",minConditionBoundary.num,";",maxConditionBoundary.num,")"))
                print(plot.gp)
            }
        # Grouped Scale(Condition & Control)
            colBreaks.num <- SuperTK::BreakVector(
                x.num=c(c(attributes(apa.mtx)$matrices$agg), c(attributes(apa.mtx)$matrices$aggCtrl)),
                n.num=51)
        # Control + Grouped Scale(with Condition)
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=attributes(apa.mtx)$matrices$aggCtrl, 
                heatmap.col=SuperTK::viridis(51),
                colBreaks.num=colBreaks.num,
                title.chr="Agregation control"
            ) + ggplot2::labs(subtitle="scale (grouped with condition)")
            print(plot.gp)
        # Condition + Grouped Scale(with Control)
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=attributes(apa.mtx)$matrices$agg, 
                heatmap.col=SuperTK::viridis(51),
                colBreaks.num=colBreaks.num,
                title.chr="Agregation"
            ) + ggplot2::labs(subtitle="scale (grouped with control)")
            print(plot.gp)
        # Grouped Scale(Corrected condition & Control)
            colBreaks.num <- SuperTK::BreakVector(
                x.num=c(c(attributes(apa.mtx)$matrices$aggCorrected), c(attributes(apa.mtx)$matrices$aggCtrl)),
                n.num=51)
        # Control + Grouped Scale(with corrected condition)
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=attributes(apa.mtx)$matrices$aggCtrl, 
                heatmap.col=SuperTK::viridis(51),
                colBreaks.num=colBreaks.num,
                title.chr="Agregation control"
            ) + ggplot2::labs(subtitle="scale (grouped with condition corrected)")
            print(plot.gp)
        # Corrected condition + Grouped Scale(with Control)
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=attributes(apa.mtx)$matrices$aggCorrected, 
                heatmap.col=SuperTK::viridis(51),
                colBreaks.num=colBreaks.num,
                title.chr="Agregation corrected"
            ) + ggplot2::labs(subtitle="scale (grouped with control)")
            print(plot.gp)
        # Density

            plot.gp <- .ggDensity(
                data.lst=list(differential=stats::na.omit(c(apa.mtx))),
                colour.col=SuperTK::Hue(3)[[1]],
                title.chr="Agregation of differential matrices density")
            print(plot.gp)

            data.lst <- list(
                control=stats::na.omit(c(attributes(apa.mtx)$matrices$aggCtrl)),
                condition=stats::na.omit(c(attributes(apa.mtx)$matrices$agg))
            )
            plot.gp <- .ggDensity(
                data.lst=data.lst,
                colour.col=SuperTK::Hue(3)[2:3],
                title.chr="Condition and control densities")
            print(plot.gp)
            
            plot.gp <- .ggDensity(
                data.lst=list(deltaCorrected=stats::na.omit(c(attributes(apa.mtx)$matrices$aggCorrectedDelta))),
                colour.col=SuperTK::Hue(3)[[1]],
                title.chr="Differential of corrected agregated matrices density")
            print(plot.gp)

            data.lst <- list(
                control=stats::na.omit(c(attributes(apa.mtx)$matrices$aggCtrl)),
                correctedCondition=stats::na.omit(c(attributes(apa.mtx)$matrices$aggCorrected))
            )
            plot.gp <- .ggDensity(
                data.lst=data.lst,
                colour.col=SuperTK::Hue(3)[2:3],
                title.chr="Condition corrected and control densities")
            print(plot.gp)
    }
    # Attributes
            grid::grid.newpage()
            attr.ndx <- apa.mtx %>%
                attributes %>%
                names %>%
                SuperTK::NotIn(c("dim","matrices","interactions", "dimnames"))
            attr.lst <- attributes(apa.mtx)[attr.ndx]
            attr.lst$aggregationMethod <- function(pxl){pxl[is.na(pxl)]<-0;mean(pxl,na.rm=TRUE)}
            attr1.ndx <- attr.lst %>%
                lapply(class) %>%
                unlist %>%
                SuperTK::NotIn(c("matrix", "list","GInteractions","function")) 
            attr1.lst <- attr.lst[attr1.ndx] %>%
                        lapply(as.character) %>%
                        unlist
            attr2.ndx <- attr.lst %>%
                        lapply(class) %>%
                        unlist %>%
                        magrittr::is_in (c("function"))
            attr2.lst <- attr.lst[attr2.ndx] %>%
                        lapply(function(function.fun){
                            function.chr <- deparse(function.fun)
                            function.chr[3:(length(function.chr)-1)] %>%
                            paste0(collapse=";\n")
                            return(gsub(" ","",function.chr))
                            }) %>%
                        unlist
            attr.lst <- c(attr1.lst, attr2.lst)
            attr.tbl  <- tibble::as_tibble_col(attr.lst) %>%
            tibble::add_column(name=names(attr.lst)) %>%
            tibble::column_to_rownames("name")
            gridExtra::grid.table(attr.tbl)
}