#' Draw multiple aggregation plot.
#'
#' PlotAPA
#' @description Draw aggregation plot from aggregation matrices.
#' @param apa.mtx <matrix>: The aggregated matrix.
#' @param trimPrct.num <numeric>: A number between 0 and 100 thaht give the percentage of triming in matrices.
#' @param colMin.num <matrix>: The minimal value in color scale. If Null automaticaly find.
#' @param colMid.num <matrix>: The middle value in color scale. If Null automaticaly find.
#' @param colMax.num <matrix>: The mximal value in color scale. If Null automaticaly find.
#' @param colCondMin.num <matrix>: Avalaible for plotting differantial aggregation. The minimal value in color scale in the classsical aggregation plot. If Null automaticaly find.
#' @param colCondMax.num <matrix>: Avalaible for plotting differantial aggregation. The maxiaml value in color scale in the classsical aggregation plot. If Null automaticaly find.
#' @return None
#' @examples
#' # Data
#' data(Beaf32_Peaks.gnr)
#' data(HiC_Ctrl.cmx_lst)
#' data(HiC_HS.cmx_lst)
#'
#' # Index Beaf32
#' Beaf32_Index.gnr <- IndexFeatures(
#'     gRange.gnr_lst = list(Beaf = Beaf32_Peaks.gnr),
#'     chromSize.dtf  = data.frame(seqnames = c("2L", "2R"), seqlengths = c(23513712, 25286936)),
#'     binSize.num    = 100000
#' )
#'
#' # Beaf32 <-> Beaf32 Pairing
#' Beaf_Beaf.gni <- SearchPairs(indexAnchor.gnr = Beaf32_Index.gnr)
#' Beaf_Beaf.gni <- Beaf_Beaf.gni[seq_len(2000)] # subset 2000 first for exemple
#'
#' # Matrices extractions center on Beaf32 <-> Beaf32 point interaction
#' interactions_Ctrl.mtx_lst <- ExtractSubmatrix(
#'     feature.gn         = Beaf_Beaf.gni,
#'     hic.cmx_lst        = HiC_Ctrl.cmx_lst,
#'     referencePoint.chr = "pf"
#' )
#' interactions_HS.mtx_lst <- ExtractSubmatrix(
#'     feature.gn         = Beaf_Beaf.gni,
#'     hic.cmx_lst        = HiC_HS.cmx_lst,
#'     referencePoint.chr = "pf"
#' )
#'
#' # Aggregate matrices in one matrix
#' aggreg.mtx <- Aggregation(interactions_Ctrl.mtx_lst)
#'
#' # MultiPlot Aggregation
#' PlotAPA(aggreg.mtx)
#'
#' # Differential Aggregation
#' aggregDiff.mtx <- Aggregation(
#'     ctrlMatrices.lst = interactions_Ctrl.mtx_lst,
#'     matrices.lst     = interactions_HS.mtx_lst,
#' )
#'
#' # MultiPlot Differential Aggregation
#' PlotAPA(aggregDiff.mtx)
#'
PlotAPA <- function(
    apa.mtx = NULL, trimPrct.num = 0, colMin.num = NULL, colMid.num = NULL,
    colMax.num = NULL, colCondMin.num = NULL, colCondMax.num = NULL
) {
    .ggDensity <- function(
        data.lst = NULL, colour.col = NULL, mean.bln = TRUE, title.chr = NULL
    ) {
        data.lst_tbl <- lapply(
            seq_along(data.lst),
            function(element.ndx) {
                return(tibble::tibble(
                    value = data.lst[[element.ndx]],
                    class = factor(names(data.lst)[[element.ndx]])
                ))
            }
        )
        data.tbl <- dplyr::bind_rows(data.lst_tbl)
        if (is.null(colour.col)) {
            colour.col <- Hue(length(data.lst)) |>
                stats::setNames(names(data.lst))
        }
        plot.gp <- ggplot2::ggplot(
            data.tbl,
            ggplot2::aes(
                x = data.tbl$value,
                fill = data.tbl$class,
                colour = data.tbl$class
            )
        ) +
            ggplot2::geom_density(alpha = 0.1) +
            ggplot2::scale_color_manual(values = colour.col) +
            ggplot2::scale_fill_manual(values = colour.col) +
            ggplot2::labs(title = title.chr)
        if (mean.bln) {
            data.tbl <- dplyr::group_by(data.tbl, class = data.tbl$class)
            mu.tbl <- dplyr::summarise(
                data.tbl,
                grp.mean = mean(data.tbl$value)
            )
            plot.gp <- plot.gp +
                ggplot2::geom_vline(
                    data = mu.tbl,
                    ggplot2::aes(
                        xintercept = mu.tbl$grp.mean,
                        colour = mu.tbl$class
                    ),
                    linetype = "dashed"
                )
        }
        return(plot.gp)
    }
    # Differential or not?
    differential.bln <- !is.null(attributes(apa.mtx)$matrices)
    if (differential.bln) {
        heatmap.col <- NULL
    } else {
        heatmap.col <- viridis(255)
    }
    # Plot Auto Scale
    plot.gp <- ggAPA(
        apa.mtx = apa.mtx,
        heatmap.col = heatmap.col,
        title.chr = ifelse(
            differential.bln,
            yes = "Agregation of differential matrices",
            no = "Agregation"
        )
    ) +
        ggplot2::labs(subtitle = "scale(auto), center()")
    plot(plot.gp)
    # Auto Scale + Center
    if (!is.null(colMid.num)) {
        plot.gp <- ggAPA(
            apa.mtx = apa.mtx,
            heatmap.col = heatmap.col,
            colMid.num = colMid.num,
            title.chr = ifelse(differential.bln,
                yes = "Agregation of differential matrices",
                no = "Agregation"
            )
        ) +
            ggplot2::labs(
                subtitle = paste0("scale(auto), center(", colMid.num, ")")
            )
        plot(plot.gp)
    }
    # Trim Scale + Center
    if (!is.null(trimPrct.num) && 0 < trimPrct.num) {
        plot.gp <- ggAPA(
            apa.mtx = apa.mtx,
            heatmap.col = heatmap.col,
            trimPrct.num = trimPrct.num,
            colMid.num = colMid.num,
            title.chr = ifelse(differential.bln,
                yes = "Agregation of differential matrices",
                no = "Agregation"
            )
        ) +
            ggplot2::labs(
                subtitle =
                    paste0(
                        "scale(rm ",trimPrct.num,"%), center(",colMid.num,")"
                    )
            )
        plot(plot.gp)
    }
    # MinMax Scale + Center
    if (!is.null(colMin.num) || !is.null(colMax.num)) {
        plot.gp <- ggAPA(
            apa.mtx = apa.mtx,
            heatmap.col = heatmap.col,
            colMin.num = colMin.num,
            colMid.num = colMid.num,
            colMax.num = colMax.num,
            title.chr = ifelse(
                differential.bln,
                yes = "Agregation of differential matrices",
                no = "Agregation"
            )
        ) +
            ggplot2::labs(
                subtitle =
                    paste0(
                        "scale(", colMin.num, ";", colMax.num, "),",
                        "center(", colMid.num, ")"
                    )
            )
        plot(plot.gp)
    }
    if (differential.bln) {
        # Pval + Auto Scale
        if (!is.null(attributes(apa.mtx)$matrices$pVal) &&
            sum(!is.na(attributes(apa.mtx)$matrices$pVal)) >= 3
        ) {
            plot.gp <- ggAPA(
                apa.mtx = attributes(apa.mtx)$matrices$pVal,
                heatmap.col = YlOrRd(9),
                title.chr = "-log10(p.values)"
            ) +
                ggplot2::labs(subtitle = "scale(auto), center()")
        } else {
            plot.gp <- ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::annotate(
                    "text",
                    x = 1,
                    y = 1,
                    label = paste0(
                        "Not enough pval computed to plot a ",
                        "pval matrix (<3) or nothing significant"
                    )
                )
        }
        plot(plot.gp)
        # FiltPval + Auto Scale + Center
        if (!is.null(attributes(apa.mtx)$matrices$aggDiffPvalFilt) &&
            sum(!is.na(attributes(apa.mtx)$matrices$pVal)) >= 3
        ) {
            plot.gp <- ggAPA(
                apa.mtx = attributes(apa.mtx)$matrices$aggDiffPvalFilt,
                heatmap.col = heatmap.col,
                colMid.num = colMid.num,
                title.chr = paste0("Agregation of differential matrices")
            ) +
                ggplot2::labs(
                    subtitle = paste0(
                        "filtred by p.values, scale(auto), center(",
                        colMid.num,")"
                    )
                )
        } else {
            plot.gp <- ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::annotate(
                    "text",
                    x = 1,
                    y = 1,
                    label = paste0(
                        "Not enough pval computed to plot a ",
                        "pval matrix (<3) or nothing significant"
                    )
                )
        }
        plot(plot.gp)
        # FiltPval + Trim Scale + Center
        if (!is.null(attributes(apa.mtx)$matrices$aggDiffPvalFilt) &&
            sum(!is.na(attributes(apa.mtx)$matrices$pVal)) >= 3
        ) {
            plot.gp <- ggAPA(
                apa.mtx = attributes(apa.mtx)$matrices$aggDiffPvalFilt,
                heatmap.col = heatmap.col,
                trimPrct.num = trimPrct.num,
                colMid.num = colMid.num,
                title.chr = paste0("Agregation of differential matrices")
            ) +
                ggplot2::labs(
                    subtitle = paste0(
                        "filtred by p.values, scale(rm ", trimPrct.num, "%), ",
                        "center(", colMid.num, ")"
                    )
                )
        } else {
            plot.gp <- ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::annotate(
                    "text",
                    x = 1,
                    y = 1,
                    label = paste0(
                        "Not enough pval computed to plot a ",
                        "pval matrix (<3) or nothing significant"
                    )
                )
        }
        plot(plot.gp)
        # Delta + Auto Scale + Center
        plot.gp <- ggAPA(
            apa.mtx = attributes(apa.mtx)$matrices$aggDelta,
            heatmap.col = heatmap.col,
            colMid.num = colMid.num,
            title.chr = "Differential of agregated matrices"
        ) +
            ggplot2::labs(
                subtitle = paste0("scale(auto), center(", colMid.num, ")")
            )
        plot(plot.gp)
        # Delta + Trim Scale + Center
        if (!is.null(trimPrct.num) && 0 < trimPrct.num) {
            plot.gp <- ggAPA(
                apa.mtx = attributes(apa.mtx)$matrices$aggDelta,
                heatmap.col = heatmap.col,
                trimPrct.num = trimPrct.num,
                colMid.num = colMid.num,
                title.chr = "Differential of agregated matrices"
            ) +
                ggplot2::labs(
                    subtitle = paste0(
                        "scale(rm ", trimPrct.num, "%), center(",colMid.num,")"
                    )
                )
            plot(plot.gp)
        }
        # Delta + Auto Scale + Center
        if (!is.null(attributes(apa.mtx)$matrices$aggCorrectedDelta)) {
            plot.gp <- ggAPA(
                apa.mtx = attributes(apa.mtx)$matrices$aggCorrectedDelta,
                heatmap.col = heatmap.col,
                colMid.num = colMid.num,
                title.chr = "Differential of corrected agregated matrices"
            ) +
                ggplot2::labs(
                    subtitle = paste0("scale (auto), center(", colMid.num, ")")
                )
            plot(plot.gp)
            # Delta + Trim Scale + Center
            if (!is.null(trimPrct.num) && 0 < trimPrct.num) {
                plot.gp <- ggAPA(
                    apa.mtx = attributes(apa.mtx)$matrices$aggCorrectedDelta,
                    heatmap.col = heatmap.col,
                    trimPrct.num = trimPrct.num,
                    colMid.num = colMid.num,
                    title.chr = "Differential of corrected agregated matrices"
                ) +
                    ggplot2::labs(
                        subtitle = paste0(
                            "scale(rm ", trimPrct.num, "%), ",
                            "center(", colMid.num, ")"
                        )
                    )
                plot(plot.gp)
            }
        }
        # Control + Auto Scale
        plot.gp <- ggAPA(
            apa.mtx = attributes(apa.mtx)$matrices$aggCtrl,
            heatmap.col = viridis(51),
            title.chr = "Agregation control"
        ) +
            ggplot2::labs(subtitle = "scale(auto)")
        plot(plot.gp)
        # Control + Trim Scale
        if (!is.null(trimPrct.num) && 0 < trimPrct.num) {
            plot.gp <- ggAPA(
                apa.mtx = attributes(apa.mtx)$matrices$aggCtrl,
                heatmap.col = viridis(51),
                trimPrct.num = trimPrct.num, title.chr = "Agregation control"
            ) +
                ggplot2::labs(subtitle = paste0("scale(rm ",trimPrct.num,"%)"))
            plot(plot.gp)
        }
        # Control + MinMax Scale
        if (!is.null(colCondMin.num) ||
            !is.null(colCondMax.num)
        ) {
            plot.gp <- ggAPA(
                apa.mtx = attributes(apa.mtx)$matrices$aggCtrl,
                heatmap.col = viridis(51),
                colMin.num = colCondMin.num, colMax.num = colCondMax.num,
                title.chr = "Agregation control"
            ) +
                ggplot2::labs(
                    subtitle = paste0(
                        "scale(", colCondMin.num, ";", colCondMax.num, ")"
                    )
                )
            plot(plot.gp)
        }
        # Condition + Auto Scale
        plot.gp <- ggAPA(
            apa.mtx = attributes(apa.mtx)$matrices$agg,
            heatmap.col = viridis(51),
            title.chr = "Agregation"
        ) +
            ggplot2::labs(subtitle = "scale(auto)")
        plot(plot.gp)
        # Condition + Trim Scale
        if (!is.null(trimPrct.num) && 0 < trimPrct.num) {
            plot.gp <- ggAPA(
                apa.mtx = attributes(apa.mtx)$matrices$agg,
                heatmap.col = viridis(51),
                trimPrct.num = trimPrct.num, title.chr = "Agregation"
            ) +
                ggplot2::labs(
                    subtitle = paste0("scale (rm ",trimPrct.num,"%)")
                )
            plot(plot.gp)
        }
        # Condition + MinMax Scale
        if (!is.null(colCondMin.num) || !is.null(colCondMax.num)) {
            plot.gp <- ggAPA(
                apa.mtx = attributes(apa.mtx)$matrices$agg,
                heatmap.col = viridis(51),
                colMin.num = colCondMin.num, colMax.num = colCondMax.num,
                title.chr = "Agregation"
            ) +
                ggplot2::labs(
                    subtitle = paste0(
                        "scale(", colCondMin.num, ";",
                        colCondMax.num, ")"
                    )
                )
            plot(plot.gp)
        }
        # Corrected condition + Auto Scale
        if (!is.null(attributes(apa.mtx)$matrices$aggCorrected)) {
            plot.gp <- ggAPA(
                apa.mtx = attributes(apa.mtx)$matrices$aggCorrected,
                heatmap.col = viridis(51),
                title.chr = "Agregation corrected"
            ) +
                ggplot2::labs(subtitle = "scale (auto)")
            plot(plot.gp)
            # Corrected condition + Trim Scale
            if (!is.null(trimPrct.num) && 0 < trimPrct.num) {
                plot.gp <- ggAPA(
                    apa.mtx = attributes(apa.mtx)$matrices$aggCorrected,
                    heatmap.col = viridis(51),
                    trimPrct.num = trimPrct.num,
                    title.chr = "Agregation corrected"
                ) +
                    ggplot2::labs(
                        subtitle = paste0("scale (rm ", trimPrct.num, "%)")
                    )
                plot(plot.gp)
            }
            # Corrected condition + MinMax Scale
            if (!is.null(colCondMin.num) || !is.null(colCondMax.num)) {
                plot.gp <- ggAPA(
                    apa.mtx = attributes(apa.mtx)$matrices$aggCorrected,
                    heatmap.col = viridis(51),
                    colMin.num = colCondMin.num,
                    colMax.num = colCondMax.num,
                    title.chr = "Agregation corrected"
                ) +
                    ggplot2::labs(
                        subtitle = paste0(
                            "scale (", colCondMin.num, ";", colCondMax.num, ")"
                        )
                    )
                plot(plot.gp)
            }
        }
        # Grouped Scale(Condition & Control)
        colBreaks.num <- BreakVector(
            x.num = c(
                c(attributes(apa.mtx)$matrices$agg),
                c(attributes(apa.mtx)$matrices$aggCtrl)
            ),
            n.num = 51
        )
        # Control + Grouped Scale(with Condition)
        plot.gp <- ggAPA(
            apa.mtx = attributes(apa.mtx)$matrices$aggCtrl,
            heatmap.col = viridis(51),
            colBreaks.num = colBreaks.num,
            title.chr = "Agregation control"
        ) +
            ggplot2::labs(subtitle = "scale (grouped with condition)")
        plot(plot.gp)
        # Condition + Grouped Scale(with Control)
        plot.gp <- ggAPA(
            apa.mtx = attributes(apa.mtx)$matrices$agg,
            heatmap.col = viridis(51),
            colBreaks.num = colBreaks.num, title.chr = "Agregation"
        ) +
            ggplot2::labs(subtitle = "scale (grouped with control)")
        plot(plot.gp)
        # Grouped Scale(Corrected condition & Control)
        if (!is.null(attributes(apa.mtx)$matrices$aggCorrected)) {
            colBreaks.num <- BreakVector(
                x.num = c(
                    c(attributes(apa.mtx)$matrices$aggCorrected),
                    c(attributes(apa.mtx)$matrices$aggCtrl)
                ),
                n.num = 51
            )
            # Control + Grouped Scale(with corrected condition)
            plot.gp <- ggAPA(
                apa.mtx = attributes(apa.mtx)$matrices$aggCtrl,
                heatmap.col = viridis(51),
                colBreaks.num = colBreaks.num,
                title.chr = "Agregation control"
            ) +
                ggplot2::labs(
                    subtitle = "scale (grouped with condition corrected)"
                )
            plot(plot.gp)
            # Corrected condition + Grouped Scale(with Control)
            plot.gp <- ggAPA(
                apa.mtx = attributes(apa.mtx)$matrices$aggCorrected,
                heatmap.col = viridis(51),
                colBreaks.num = colBreaks.num,
                title.chr = "Agregation corrected"
            ) +
                ggplot2::labs(subtitle = "scale (grouped with control)")
            plot(plot.gp)
        }
        # Density
        plot.gp <- .ggDensity(
            data.lst = list(differential = stats::na.omit(c(apa.mtx))),
            colour.col = Hue(3)[[1]],
            title.chr = "Agregation of differential matrices density"
        )
        plot(plot.gp)

        data.lst <- list(
            control = stats::na.omit(c(attributes(apa.mtx)$matrices$aggCtrl)),
            condition = stats::na.omit(c(attributes(apa.mtx)$matrices$agg))
        )
        plot.gp <- .ggDensity(
            data.lst = data.lst, colour.col = Hue(3)[2:3],
            title.chr = "Condition and control densities"
        )
        plot(plot.gp)
        if (!is.null(attributes(apa.mtx)$matrices$aggCorrectedDelta)) {
            plot.gp <- .ggDensity(
                data.lst = list(
                    deltaCorrected = stats::na.omit(
                        c(attributes(apa.mtx)$matrices$aggCorrectedDelta)
                    )
                ),
                colour.col = Hue(3)[[1]],
                title.chr = paste0(
                    "Differential of corrected ",
                    "agregated matrices density"
                )
            )
            plot(plot.gp)
        }
        if (!is.null(attributes(apa.mtx)$matrices$aggCorrected)) {
            data.lst <- list(
                control = stats::na.omit(
                    c(attributes(apa.mtx)$matrices$aggCtrl)
                ),
                correctedCondition = stats::na.omit(
                    c(attributes(apa.mtx)$matrices$aggCorrected)
                )
            )
            plot.gp <- .ggDensity(
                data.lst = data.lst, colour.col = Hue(3)[2:3],
                title.chr = "Condition corrected and control densities"
            )
            plot(plot.gp)
        }
    }
    # Attributes
    grid::grid.newpage()
    attr.ndx <- apa.mtx |>
        attributes() |>
        names() |>
        NotIn(
            c("dim", "matrices", "interactions", "dimnames", "correctionArea")
        )
    attr.lst <- attributes(apa.mtx)[attr.ndx]
    attr1.ndx <- attr.lst |>
        lapply(class) |>
        unlist() |>
        NotIn(c("matrix", "list", "GInteractions", "function"))
    attr1.lst <- attr.lst[attr1.ndx] |>
        lapply(
            function(attr) {
                as.character(attr) |>
                    paste0(collapse = ",")
            }
        ) |>
        unlist()
    attr2.ndx <- unlist(lapply(attr.lst, class)) %in% "function"
    attr2.lst <- attr.lst[attr2.ndx] |>
        lapply(
            function(fun) {
                fun.chr <- deparse(fun)
                fun.chr <- fun.chr[3:(length(fun.chr) - 1)] |>
                    paste0(collapse = ";\n")
                return(gsub(" ", "", fun.chr))
            }
        ) |>
        unlist()
    attr.lst <- c(attr1.lst, attr2.lst)
    attr.tbl <- tibble::as_tibble_col(attr.lst) |>
        tibble::add_column(name = names(attr.lst)) |>
        tibble::column_to_rownames("name")
    gridExtra::grid.table(attr.tbl)
}
