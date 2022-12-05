#' Compute HiC matrix-balancing.
#'
#' NormalizeHiC
#' @description Apply a matrix-balancing normalization method to a list of contacts matrix.
#' @param hic.cmx_lst <List[contactMatrix]>: The HiC maps list.
#' @param method.chr <character> : The kind of normalization method. One of "ICE", "VC" or "VC_SQRT" (Default "ICE")
#' @param interaction.type <character> : "cis", "trans", c("cis", "trans"), "all". If NULL normalization is apply on cis chunks then trans chunks (equivalent to c("cis", "trans")). If is "all", normalization is apply on all chuncks at once. (Default NULL)
#' @param maxIter.num <numerical>: The maximum iteration number.
#' @param qtlTh.num <numerical>: the threshold quantile below which the bins will be ignored. (Default 0.15)
#' @param cores.num <numerical> : Number of cores to be used. (Default 1)
#' @param verbose.bln <logical>: If TRUE show the progression in console. (Default FALSE)
#' @return A matrices list.
#' @examples
#' data(HiC_Ctrl.cmx_lst)
#'
#' HiC_Ctrl_ICE.cmx_lst <- NormalizeHiC(HiC_Ctrl.cmx_lst,
#'     interaction.type = "cis",
#'     method.chr = "ICE"
#' )
#'
#' HiC_Ctrl_VC.cmx_lst <- NormalizeHiC(HiC_Ctrl.cmx_lst,
#'     interaction.type = c("cis", "trans"),
#'     method.chr = "VC"
#' )
#'
#' HiC_Ctrl_VC_SQRT.cmx_lst <- NormalizeHiC(HiC_Ctrl.cmx_lst,
#'     interaction.type = "all",
#'     method.chr = "VC_SQRT"
#' )
#'
NormalizeHiC <- function(
	hic.cmx_lst, method.chr = "ICE", interaction.type = NULL, maxIter.num = 50,
	qtlTh.num = 0.15, cores.num = 1, verbose.bln = FALSE
) {
	if (!is.null(interaction.type) &&
		"all" %in% interaction.type
	) {
		megaHic.cmx <- JoinHiC(hic.cmx_lst)
		if (method.chr == "VC") {
			megaHic.cmx <- VCnorm(
				megaHic.cmx,
				qtlTh.num = qtlTh.num,
				sqrt.bln = FALSE
			)
		} else if (method.chr == "VC_SQRT") {
			megaHic.cmx <- VCnorm(
				megaHic.cmx,
				qtlTh.num = qtlTh.num,
				sqrt.bln = TRUE
			)
		} else if (method.chr == "ICE") {
			megaHic.cmx <- ICEnorm(
				megaHic.cmx,
				qtlTh.num = qtlTh.num,
				maxIter.num = maxIter.num
			)
		}
		hic.cmx_lst <- CutHiC(megaHic.cmx, verbose.bln = verbose.bln)
	} else if (!is.null(interaction.type) &&
		"cis" %in% interaction.type &&
		NotIn("trans", interaction.type)
	) {
		matricesKind.tbl <- attributes(hic.cmx_lst)$matricesKind
		if ("cis" %in% matricesKind.tbl$type) {
			cisMatricesNames.chr <- dplyr::filter(
					matricesKind.tbl,
					matricesKind.tbl$type == "cis"
				) |>
				dplyr::pull("name")
			if (length(cisMatricesNames.chr)) {
				multicoreParam <- MakeParallelParam(
					cores.num = cores.num,
					verbose.bln = verbose.bln
				)
				hic.cmx_lst[cisMatricesNames.chr] <- BiocParallel::bplapply(
					BPPARAM = multicoreParam, seq_along(cisMatricesNames.chr),
					function(ele.ndx) {
						matrixName.chr <- cisMatricesNames.chr[[ele.ndx]]
						if (method.chr == "VC") {
							hic.cmx <- VCnorm(
								hic.cmx_lst[[matrixName.chr]],
								qtlTh.num = qtlTh.num,
								sqrt.bln = FALSE
							)
						} else if (method.chr == "VC_SQRT") {
							hic.cmx <- VCnorm(
								hic.cmx_lst[[matrixName.chr]],
								qtlTh.num = qtlTh.num,
								sqrt.bln = TRUE
							)
						} else if (method.chr == "ICE") {
							hic.cmx <- ICEnorm(
								hic.cmx_lst[[matrixName.chr]],
								qtlTh.num = qtlTh.num,
								maxIter.num = maxIter.num
							)
						}
						return(hic.cmx)
					}
				)
			}
			transMatricesNames.chr <- dplyr::filter(
				matricesKind.tbl,
				matricesKind.tbl$type == "trans") |>
				dplyr::pull("name")
			message(paste0(
				paste(transMatricesNames.chr, collapse = ", "),
				" remove from output."
			))
			if (length(transMatricesNames.chr)) {
				attr.lst <- attributes(hic.cmx_lst)
				attr.lst$matricesKind <- dplyr::filter(
					attr.lst$matricesKind,
					NotIn(attr.lst$matricesKind$name, transMatricesNames.chr)
				)
				chroms.chr <- attr.lst$matricesKind$name |>
					strsplit("_") |>
					unlist() |>
					unique()
				attr.lst$chromSize <- dplyr::filter(
					attr.lst$chromSize,
					attr.lst$chromSize$name == chroms.chr
				)
				hic.cmx_lst <- hic.cmx_lst[-which(
					names(hic.cmx_lst) %in% transMatricesNames.chr
				)] |>
					AddAttr(attr.lst)
			}
		} else {
			message(paste0(
				"No cis matrix, Normalisation ",
				"won't be applied on cis matrices"
			))
		}
	} else if (!is.null(interaction.type) &&
		"trans" %in% interaction.type && NotIn("cis", interaction.type)
	) {
		matricesKind.tbl <- attributes(hic.cmx_lst)$matricesKind
		if ("trans" %in% matricesKind.tbl$type) {
			transMatricesNames.chr <- dplyr::filter(
					matricesKind.tbl,
					matricesKind.tbl$type == "trans") |>
				dplyr::pull("name")
			chromNames.chr <- transMatricesNames.chr |>
				strsplit("_") |>
				unlist() |>
				unique()
			chromSize.tbl <- attributes(hic.cmx_lst)$chromSize
			chromSize.tbl <- dplyr::filter(
				chromSize.tbl,
				chromSize.tbl$name %in% chromNames.chr
			)
			trans.cmx_lst <- hic.cmx_lst[transMatricesNames.chr] |>
				AddAttr(list(
					resolution = attributes(hic.cmx_lst)$resolution,
					chromSize = chromSize.tbl,
					matricesKind = matricesKind.tbl
				))
			megaHic.cmx <- JoinHiC(trans.cmx_lst)
			if (method.chr == "VC") {
				megaHic.cmx <- VCnorm(
					megaHic.cmx,
					qtlTh.num = qtlTh.num,
					sqrt.bln = FALSE
				)
			} else if (method.chr == "VC_SQRT") {
				megaHic.cmx <- VCnorm(
					megaHic.cmx,
					qtlTh.num = qtlTh.num,
					sqrt.bln = TRUE
				)
			} else if (method.chr == "ICE") {
				megaHic.cmx <- ICEnorm(
					megaHic.cmx,
					qtlTh.num = qtlTh.num,
					maxIter.num = maxIter.num
				)
			}
			trans.cmx_lst <- CutHiC(megaHic.cmx, verbose.bln = verbose.bln)
			hic.cmx_lst[transMatricesNames.chr] <-
				trans.cmx_lst[transMatricesNames.chr]
			cisMatricesNames.chr <- dplyr::filter(
					matricesKind.tbl,
					matricesKind.tbl$type == "cis"
				) |>
				dplyr::pull("name")
			message(paste0(
				paste(cisMatricesNames.chr, collapse = ", "),
				" remove from output."
			))
			if (length(cisMatricesNames.chr)) {
				attr.lst <- attributes(hic.cmx_lst)
				attr.lst$matricesKind <- dplyr::filter(
					attr.lst$matricesKind,
					NotIn(attr.lst$matricesKind$name, cisMatricesNames.chr)
			 	)
				chroms.chr <- attr.lst$matricesKind$name |>
					strsplit("_") |>
					unlist() |>
					unique()
				attr.lst$chromSize <- dplyr::filter(
					attr.lst$chromSize,
					attr.lst$chromSize$name == chroms.chr
				)
				hic.cmx_lst <- hic.cmx_lst[-which(
					names(hic.cmx_lst) %in% cisMatricesNames.chr
				)] |>
					AddAttr(attr.lst)
			}
		} else {
			message(paste0(
				"No trans matrix, Normalisation ",
				"won't be applied on trans matrices"
			))
		}
	} else {
		matricesKind.tbl <- attributes(hic.cmx_lst)$matricesKind
		if (is.null(interaction.type) |
			"cis" %in% interaction.type
		) {
			if ("cis" %in% matricesKind.tbl$type) {
				cisMatricesNames.chr <- dplyr::filter(
						matricesKind.tbl,
						matricesKind.tbl$type == "cis"
					) |>
					dplyr::pull("name")
				if (length(cisMatricesNames.chr)) {
					multicoreParam <- MakeParallelParam(
						cores.num = cores.num,
						verbose.bln = verbose.bln
					)
					hic.cmx_lst[cisMatricesNames.chr] <-
						BiocParallel::bplapply(
							BPPARAM = multicoreParam,
							seq_along(cisMatricesNames.chr),
							function(ele.ndx) {
							matrixName.chr <- cisMatricesNames.chr[[ele.ndx]]
							if (method.chr == "VC") {
								hic.cmx <- VCnorm(
									hic.cmx_lst[[matrixName.chr]],
									qtlTh.num = qtlTh.num,
									sqrt.bln = FALSE
								)
							} else if (method.chr == "VC_SQRT") {
								hic.cmx <- VCnorm(
									hic.cmx_lst[[matrixName.chr]],
									qtlTh.num = qtlTh.num,
									sqrt.bln = TRUE
								)
							} else if (method.chr == "ICE") {
								hic.cmx <- ICEnorm(
									hic.cmx_lst[[matrixName.chr]],
									qtlTh.num = qtlTh.num,
									maxIter.num = maxIter.num
								)
							}
							return(hic.cmx)
							}
					)
				}
			} else {
				message(paste0(
					"No cis matrix, Normalisation ",
					"won't be applied on cis matrices"
				))
			}
		}
		if (is.null(interaction.type) |
			"trans" %in% interaction.type
		) {
			if ("trans" %in% matricesKind.tbl$type) {
				transMatricesNames.chr <- dplyr::filter(
						matricesKind.tbl,
						matricesKind.tbl$type == "trans"
					) |>
					dplyr::pull("name")
				chromNames.chr <- transMatricesNames.chr |>
					strsplit("_") |>
					unlist() |>
					unique()
				chromSize.tbl <- attributes(hic.cmx_lst)$chromSize
				chromSize.tbl <- dplyr::filter(
					chromSize.tbl,
					chromSize.tbl$name %in% chromNames.chr
				)
				trans.cmx_lst <- hic.cmx_lst[transMatricesNames.chr] |>
					AddAttr(list(
						resolution = attributes(hic.cmx_lst)$resolution,
						chromSize = chromSize.tbl,
						matricesKind = matricesKind.tbl
					))
				megaHic.cmx <- JoinHiC(trans.cmx_lst)
				if (method.chr == "VC") {
					megaHic.cmx <- VCnorm(
						megaHic.cmx,
						qtlTh.num = qtlTh.num,
						sqrt.bln = FALSE
					)
				} else if (method.chr == "VC_SQRT") {
					megaHic.cmx <- VCnorm(
						megaHic.cmx,
						qtlTh.num = qtlTh.num,
						sqrt.bln = TRUE
					)
				} else if (method.chr == "ICE") {
					megaHic.cmx <- ICEnorm(
						megaHic.cmx,
						qtlTh.num = qtlTh.num,
						maxIter.num = maxIter.num
					)
				}
				trans.cmx_lst <- CutHiC(megaHic.cmx, verbose.bln = verbose.bln)
				hic.cmx_lst[transMatricesNames.chr] <-
					trans.cmx_lst[transMatricesNames.chr]
			} else {
			message(paste0(
				"No trans matrix, Normalisation ",
				"won't be applied on trans matrices"
			))
			}
		}
	}
	return(hic.cmx_lst)
}
