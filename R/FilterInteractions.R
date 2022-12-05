#' Submatrix or Interactions filtering.
#'
#' FilterInteractions
#' @description Search in a GInteraction object which interactions correspond ti a target list and return a list of index or filter a matrices list according to target and a selection function.
#' @param matrices.lst <List[matrix]>: The matrices list to filter. If is not NULL, the function will return the filtred matrices list, else return a list of index.
#' @param interarctions.gni <GInteractions>: The GInteraction object on which compute the filter.
#' @param target.lst <List>: a nammed list that describe the target.
#' @param selection.fun <function>: A function that described how the target must be cross. (Defaul intersection of all targets)
#' @return A list of elements index or a filtred matrices list with attributes updates.
#' @examples
#' # Data
#' data(Beaf32_Peaks.gnr)
#' data(HiC_Ctrl.cmx_lst)
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
#' interactions_PF.mtx_lst  <- ExtractSubmatrix(
#'     feature.gn         = Beaf_Beaf.gni,
#'     hic.cmx_lst        = HiC_Ctrl.cmx_lst,
#'     referencePoint.chr = "pf"
#' )
#'
#' # Create a target
#' target.lst <- list(
#'     anchor.Beaf.name = c("Beaf32_108", "Beaf32_814"),
#'     distance         = function(dist) {
#'         dist < 300000
#'     }
#' )
#' # We target the Beaf32<->Beaf32 interactions that are less than 300Kb away
#' # and have peak Beaf32_2 and Beaf32_191 as anchors (i.e. left).
#'
#' # Create a selection
#' selection.fun <- function() {
#'     intersect(anchor.Beaf.name, distance)
#' }
#' # We select the Beaf32<->Beaf32 interactions that satisfy both targeting
#' # criteria (intersection).
#'
#' # Filtration on InteractionSet (Beaf32 <-> Beaf32 Pairs)
#' FilterInteractions(
#'     interarctions.gni = Beaf_Beaf.gni,
#'     target.lst        = target.lst,
#'     selection.fun     = NULL
#' ) |> str(max.level = 1)
#' # Returns a named list (the names match the targeting criteria).
#' # Each element is an index vector of Beaf32<->Beaf32 interactions
#' # that satisfy the given criteria.
#'
#' # Filtration on Matrices List (Beaf32 <-> Beaf32 Extracted matrices)
#' FilterInteractions(
#'     matrices.lst      = interactions_PF.mtx_lst,
#'     target.lst        = target.lst,
#'     selection.fun     = NULL
#' ) |> str(max.level = 1)
#' # Return the same kind of result.
#'
#' # Add the selection on InteractionSet Filtration
#' FilterInteractions(
#'     interarctions.gni = Beaf_Beaf.gni,
#'     target.lst        = target.lst,
#'     selection.fun     = selection.fun
#' ) |> str(max.level = 1)
#' # This return the intersection of the index vector that satisfy both
#' # targeting criteria.
#'
#' # Add the selection on Matrices List Filtration
#' FilterInteractions(
#'     matrices.lst      = interactions_PF.mtx_lst,
#'     target.lst        = target.lst,
#'     selection.fun     = selection.fun
#' ) |> str(max.level = 1)
#' # This return the filtred matrices.lst, i.e the matrices.lst for which
#' # the Beaf32<->Beaf32 interactions satisfy both targeting criteria.
#'
#' # Filtration with InteractionsSet as filtration criteria
#' target.lst <- list(interactions = Beaf_Beaf.gni[seq_len(2)])
#' FilterInteractions(
#'     interarctions.gni = Beaf_Beaf.gni,
#'     target.lst        = target.lst,
#'     selection.fun     = NULL
#' ) |> str(max.level = 1)
#'
#'
#' # Filtration with GRanges as filtration criteria
#' target.lst <- list(first =  InteractionSet::anchors(Beaf_Beaf.gni)[["first"]][seq_len(2)])
#' FilterInteractions(
#'     interarctions.gni = Beaf_Beaf.gni,
#'     target.lst        = target.lst,
#'     selection.fun     = NULL
#' ) |> str(max.level = 1)
#'
FilterInteractions <- function(
    matrices.lst = NULL, interarctions.gni = NULL, target.lst = NULL,
    selection.fun = function() {Reduce(intersect, interarctions.ndx_lst)}
) {
    if (!is.null(matrices.lst) &&
        !is.null(attributes(matrices.lst)$interactions)) {
        interarctions.gni <- attributes(matrices.lst)$interactions
    }
    interarctions.ndx_lst <- lapply(
        seq_along(target.lst),
        function(target.ndx) {
            columnName.chr <- names(target.lst)[target.ndx]
            if (is.function(target.lst[[target.ndx]])) {
                lapply(
                    S4Vectors::mcols(interarctions.gni)[,columnName.chr],
                    target.lst[[target.ndx]]
                ) |>
                    unlist() |>
                    which()
            } else if (is.character(target.lst[[target.ndx]])) {
                lapply(
                    S4Vectors::mcols(interarctions.gni)[,columnName.chr],
                    function(columnElement) {
                        res.lgk <- intersect(
                        as.character(columnElement),
                        target.lst[[columnName.chr]]
                        ) |>
                        length() |>
                        as.logical()
                        return(res.lgk)
                    }
                ) |>
                    unlist() |>
                    which()
            } else if (inherits(target.lst[[target.ndx]], "GRanges")) {
                GenomicRanges::findOverlaps(
                    InteractionSet::anchors(
                        interarctions.gni)[[columnName.chr]],
                    target.lst[[target.ndx]]
                )@from
            } else if (inherits(target.lst[[target.ndx]], "GInteractions")) {
                InteractionSet::findOverlaps(
                    interarctions.gni,
                    target.lst[[target.ndx]]
                )@from
            }
        }
    ) |>
        stats::setNames(names(target.lst))
    if (length(target.lst) == 1) {
        interarctions.ndx <- unlist(interarctions.ndx_lst)
    } else if (!is.null(selection.fun)) {
        for (target.ndx in seq_along(interarctions.ndx_lst)) {
            assign(
                names(interarctions.ndx_lst)[target.ndx],
                interarctions.ndx_lst[[target.ndx]],
                envir = parent.frame()
            )
        }
        interarctions.ndx <- selection.fun()
    } else {
        return(interarctions.ndx_lst)
    }
    if (!is.null(matrices.lst)) {
        matrices.filt.lst <- matrices.lst[interarctions.ndx]
        attributes(matrices.filt.lst)$interactions <-
            attributes(matrices.lst)$interactions[interarctions.ndx]
        attributes(matrices.filt.lst)$resolution <-
            attributes(matrices.lst)$resolution
        attributes(matrices.filt.lst)$referencePoint <-
            attributes(matrices.lst)$referencePoint
        attributes(matrices.filt.lst)$matriceDim <-
            attributes(matrices.lst)$matriceDim
        attributes(matrices.filt.lst)$target <- target.lst
        attributes(matrices.filt.lst)$selection <- selection.fun
        return(matrices.filt.lst)
    } else {
        return(interarctions.ndx)
    }
}
