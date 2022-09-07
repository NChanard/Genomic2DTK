#' DrawAPA
#'
#' Draw aggregation plot from aggregation matrices.
#' @param apa.mtx <matrix>: The aggregated matrix.
#' @param fileName.chr <character>: A string to name the ouput pdf file. If NULL is automatically generated with Versioning() function from DevTK package.
#' @param out.dir <character>: The directory where the output file will be save. (Default working directory)
#' @param trimPrct.num <numeric>: A number between 0 and 100 thaht give the percentage of triming in matrices.
#' @param minBoundary.num <matrix>: The minimal value in color scale. If Null automaticaly find.
#' @param center.num <matrix>: The middle value in color scale. If Null automaticaly find.
#' @param maxBoundary.num <matrix>: The mximal value in color scale. If Null automaticaly find.
#' @param minConditionBoundary.num <matrix>: Avalaible for plotting differantial aggregation. The minimal value in color scale in the classsical aggregation plot. If Null automaticaly find.
#' @param maxConditionBoundary.num <matrix>: Avalaible for plotting differantial aggregation. The maxiaml value in color scale in the classsical aggregation plot. If Null automaticaly find.
#' @return None
DrawAPA = function(apa.mtx = NULL, fileName.chr = DevTK::Versioning(), out.dir = getwd(), trimPrct.num=0, minBoundary.num=NULL, center.num=NULL, maxBoundary.num=NULL, minConditionBoundary.num=NULL, maxConditionBoundary.num=NULL){
    .ggDensity <- function(data.lst=NULL, colour.col=NULL, mean.bln=TRUE, title.chr=NULL){
        data.tbl <- lapply(seq_along(data.lst),function(element.ndx){
            tibble::tibble(
                value = data.lst[[element.ndx]],
                class = factor(names(data.lst)[[element.ndx]])
            ) %>%
            return(.)
            }) %>%
        dplyr::bind_rows(.)
        if(is.null(colour.col)){
            colour.col <- ColorTK::Hue(length(data.lst)) %>% magrittr::set_names(names(data.lst))
        }
        plot.gp <- ggplot2::ggplot(data.tbl, ggplot2::aes(x=value, fill=class, colour=class)) + 
            ggplot2::geom_density(alpha=0.1) +
            ggplot2::scale_color_manual(values = colour.col)+
            ggplot2::scale_fill_manual(values = colour.col) +
            ggplot2::labs(title=title.chr) 
        if (mean.bln){
            mu.tbl <- data.tbl %>% dplyr::group_by(class) %>% dplyr::summarise(grp.mean = mean(value))
            plot.gp <- plot.gp + 
                ggplot2::geom_vline(data = mu.tbl, ggplot2::aes(xintercept = grp.mean, colour = class), linetype = "dashed")
        }
        return(plot.gp)
    }
    # Differential or not?
        differential.bln <- !is.null(attributes(apa.mtx)$matrices)
    # File Path
        fileName.chr <- differential.bln %>%
            ifelse(., yes="APA_diff_", no="APA_") %>%
            paste0(.,fileName.chr)
        if (dir.exists(out.dir) && (stringr::str_ends(out.dir, "/", negate=TRUE))) {out.dir %<>% paste0("/")}
        file.pth <- paste0(out.dir, fileName.chr,".pdf")
    # Colors
        if(differential.bln){
            heatmap.col = NULL
        }else{
            heatmap.col = viridisLite::viridis(11)
        }
    # Plot
        grDevices::pdf(file.pth)
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
                    heatmap.col=RColorBrewer::brewer.pal(9, 'YlOrRd'),
                    title.chr = "-log10(p.values)"
                ) + ggplot2::labs(subtitle="scale (auto), center()")
            }else{
                plot.gp <- ggplot2::ggplot() +
                    ggplot2::theme_void() +
                    ggplot2::annotate("text", x = 1, y = 1, label = "Not enough pval computed to plot a pval matrix (<3) or nothing significant")
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
                    ggplot2::annotate("text", x = 1, y = 1, label = "Not enough pval computed to plot a pval matrix (<3) or nothing significant")
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
                    ggplot2::annotate("text", x = 1, y = 1, label = "Not enough pval computed to plot a pval matrix (<3) or nothing significant")
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
                heatmap.col=viridisLite::viridis(51),
                title.chr="Agregation control"
            ) + ggplot2::labs(subtitle="scale (auto)")
            print(plot.gp)
        # Control + Trim Scale
            if(!is.null(trimPrct.num) && 0<trimPrct.num){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$aggCtrl, 
                    heatmap.col=viridisLite::viridis(51),
                    trimPrct.num=trimPrct.num,
                    title.chr="Agregation control"
                ) + ggplot2::labs(subtitle=paste0("scale (rm ",trimPrct.num,"%)"))
                print(plot.gp)
            }
        # Control + MinMax Scale
            if(!is.null(minConditionBoundary.num) || !is.null(maxConditionBoundary.num)){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$aggCtrl, 
                    heatmap.col=viridisLite::viridis(51),
                    minBoundary.num=minConditionBoundary.num,
                    maxBoundary.num=maxConditionBoundary.num,
                    title.chr="Agregation control"
                ) + ggplot2::labs(subtitle=paste0("scale (",minConditionBoundary.num,";",maxConditionBoundary.num,")"))
                print(plot.gp)
            }
        # Condition + Auto Scale
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=attributes(apa.mtx)$matrices$agg, 
                heatmap.col=viridisLite::viridis(51),
                title.chr="Agregation"
            ) + ggplot2::labs(subtitle="scale (auto)")
            print(plot.gp)
        # Condition + Trim Scale
            if(!is.null(trimPrct.num) && 0<trimPrct.num){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$agg, 
                    heatmap.col=viridisLite::viridis(51),
                    trimPrct.num=trimPrct.num,
                    title.chr="Agregation"
                ) + ggplot2::labs(subtitle=paste0("scale (rm ",trimPrct.num,"%)"))
                print(plot.gp)
            }
        # Condition + MinMax Scale
            if(!is.null(minConditionBoundary.num) || !is.null(maxConditionBoundary.num)){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$agg, 
                    heatmap.col=viridisLite::viridis(51),
                    minBoundary.num=minConditionBoundary.num,
                    maxBoundary.num=maxConditionBoundary.num,
                    title.chr="Agregation"
                ) + ggplot2::labs(subtitle=paste0("scale (",minConditionBoundary.num,";",maxConditionBoundary.num,")"))
                print(plot.gp)
            }
        # Corrected condition + Auto Scale
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=attributes(apa.mtx)$matrices$aggCorrected, 
                heatmap.col=viridisLite::viridis(51),
                title.chr="Agregation corrected"
            ) + ggplot2::labs(subtitle="scale (auto)")
            print(plot.gp)
        # Corrected condition + Trim Scale
            if(!is.null(trimPrct.num) && 0<trimPrct.num){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$aggCorrected, 
                    heatmap.col=viridisLite::viridis(51),
                    trimPrct.num=trimPrct.num,
                    title.chr="Agregation corrected"
                ) + ggplot2::labs(subtitle=paste0("scale (rm ",trimPrct.num,"%)"))
                print(plot.gp)
            }
        # Corrected condition + MinMax Scale
            if(!is.null(minConditionBoundary.num) || !is.null(maxConditionBoundary.num)){
                plot.gp <- Genomic2DTK::ggAPA(
                    apa.mtx=attributes(apa.mtx)$matrices$aggCorrected, 
                    heatmap.col=viridisLite::viridis(51),
                    minBoundary.num=minConditionBoundary.num,
                    maxBoundary.num=maxConditionBoundary.num,
                    title.chr="Agregation corrected"
                ) + ggplot2::labs(subtitle=paste0("scale (",minConditionBoundary.num,";",maxConditionBoundary.num,")"))
                print(plot.gp)
            }
        # Grouped Scale(Condition & Control)
            colBreaks.num <- c(c(attributes(apa.mtx)$matrices$agg), c(attributes(apa.mtx)$matrices$aggCtrl)) %>% StatTK::BreakVector(x.num=., n.num=51)
        # Control + Grouped Scale(with Condition)
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=attributes(apa.mtx)$matrices$aggCtrl, 
                heatmap.col=viridisLite::viridis(51),
                colBreaks.num=colBreaks.num,
                title.chr="Agregation control"
            ) + ggplot2::labs(subtitle="scale (grouped with condition)")
            print(plot.gp)
        # Condition + Grouped Scale(with Control)
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=attributes(apa.mtx)$matrices$agg, 
                heatmap.col=viridisLite::viridis(51),
                colBreaks.num=colBreaks.num,
                title.chr="Agregation"
            ) + ggplot2::labs(subtitle="scale (grouped with control)")
            print(plot.gp)
        # Grouped Scale(Corrected condition & Control)
            colBreaks.num <- c(c(attributes(apa.mtx)$matrices$aggCorrected), c(attributes(apa.mtx)$matrices$aggCtrl)) %>% StatTK::BreakVector(x.num=., n.num=51)
        # Control + Grouped Scale(with corrected condition)
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=attributes(apa.mtx)$matrices$aggCtrl, 
                heatmap.col=viridisLite::viridis(51),
                colBreaks.num=colBreaks.num,
                title.chr="Agregation control"
            ) + ggplot2::labs(subtitle="scale (grouped with condition corrected)")
            print(plot.gp)
        # Corrected condition + Grouped Scale(with Control)
            plot.gp <- Genomic2DTK::ggAPA(
                apa.mtx=attributes(apa.mtx)$matrices$aggCorrected, 
                heatmap.col=viridisLite::viridis(51),
                colBreaks.num=colBreaks.num,
                title.chr="Agregation corrected"
            ) + ggplot2::labs(subtitle="scale (grouped with control)")
            print(plot.gp)
        # Density

            plot.gp <- .ggDensity(data.lst=list(differential=stats::na.omit(c(apa.mtx))), colour.col=ColorTK::Hue(3)[[1]], title.chr="Agregation of differential matrices density")
            print(plot.gp)

            data.lst <- list(
                control=stats::na.omit(c(attributes(apa.mtx)$matrices$aggCtrl)),
                condition=stats::na.omit(c(attributes(apa.mtx)$matrices$agg))
            )
            plot.gp <- .ggDensity(data.lst=data.lst, colour.col=ColorTK::Hue(3)[2:3], title.chr="Condition and control densities")
            print(plot.gp)
            
            plot.gp <- .ggDensity(data.lst=list(deltaCorrected=stats::na.omit(c(attributes(apa.mtx)$matrices$aggCorrectedDelta))), colour.col=ColorTK::Hue(3)[[1]], title.chr="Differential of corrected agregated matrices density")
            print(plot.gp)

            data.lst <- list(
                control=stats::na.omit(c(attributes(apa.mtx)$matrices$aggCtrl)),
                correctedCondition=stats::na.omit(c(attributes(apa.mtx)$matrices$aggDelta))
            )
            plot.gp <- .ggDensity(data.lst=data.lst, colour.col=ColorTK::Hue(3)[2:3], title.chr="Condition corrected and control densities")
            print(plot.gp)
    }
    # Attributes
            grid::grid.newpage()
            attr.lst <- apa.mtx %>% attributes %>% names %>% DevTK::NotIn(.,c("dim","matrices","interactions", "dimnames")) %>% attributes(apa.mtx)[.]
            attr.lst$aggregationMethod <- function(pxl){pxl[is.na(pxl)]<-0;mean(pxl,na.rm=TRUE)}

            attr.lst <- c(
                    {attr.lst %>% lapply(., class) %>% unlist %>% DevTK::NotIn(.,c("matrix", "list","GInteractions","function")) %>% attr.lst[.] %>% lapply(.,as.character) %>% unlist},
                    {attr.lst %>% lapply(., class) %>% unlist %>%  magrittr::is_in (.,c("function")) %>% attr.lst[.] %>% lapply(.,function(function.fun){function.fun %>% deparse(.) %>%.[3:(length(.)-1)] %>% paste0(.,collapse=";\n") %>% stringr::str_remove_all(.," ") %>% return(.)}) %>% unlist}
                )
            attr.tbl  <- tibble::as_tibble_col(attr.lst) %>% tibble::add_column(name=names(attr.lst)) %>%  tibble::column_to_rownames(.,"name")
            gridExtra::grid.table(attr.tbl)
            grDevices::dev.off()
}