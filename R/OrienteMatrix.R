#' Extracted matrix orientation
#'
#' OrienteMatrix
#' @description Oriente extracted Matrix according to the anchors and bait order. Apply a 180° rotation follow with a transposation on a matrix or on matricies in a list according to the interactions attributes of the list.
#' @param matrice.mtx <matrix or List[matrix]>: matrix or matricies list to oriente
#' @return Oriented matrix or matricies list
#' @examples
#' # Data
#' data(Beaf32_Peaks.gnr)
#' data(HiC_Ctrl.cmx_lst)
#'
#' # Index Beaf32 in TADs domains
#' Beaf32_Index.gnr <- IndexFeatures(
#'     gRange.gnr_lst = list(Beaf = Beaf32_Peaks.gnr),
#'     chromSize.dtf = data.frame(
#'         seqnames = c("2L", "2R"),
#'         seqlengths = c(23513712, 25286936)
#'     ),
#'     binSize.num = 100000
#' )
#'
#' # Beaf32 <-> Beaf32 Pairing
#' Beaf_Beaf.gni <- SearchPairs(indexAnchor.gnr = Beaf32_Index.gnr)
#' Beaf_Beaf.gni <- Beaf_Beaf.gni[seq_len(2000)] # subset 2000 first for exemple
#'
#' # Matrices extractions center on Beaf32 <-> Beaf32 point interaction
#' interactions_PF.mtx_lst <- ExtractSubmatrix(
#'     feature.gn         = Beaf_Beaf.gni,
#'     hic.cmx_lst        = HiC_Ctrl.cmx_lst,
#'     referencePoint.chr = "pf"
#' )
#'
#' # Matrices Orientation
#' oriented_Interactions_PF.mtx_lst <- OrienteMatrix(interactions_PF.mtx_lst)
#'
OrienteMatrix <- function(
    matrice.mtx
) {
    if (is.list(matrice.mtx) &&
        !is.null(attributes(matrice.mtx)$interactions)
    ) {
        orientedMatrice.mtx <- matrice.mtx
        message(
            length(which(!attributes(matrice.mtx)$interactions$orientation)),
            " matrices are oriented"
        )
        orientedMatrice.mtx[which(
            !attributes(matrice.mtx)$interactions$orientation
        )] <- lapply(
            orientedMatrice.mtx[which(
                !attributes(matrice.mtx)$interactions$orientation
            )],
            OrienteMatrix
        )
        orientedMatrice.mtx <- AddAttr(
            orientedMatrice.mtx,
            attributes(matrice.mtx)
        )
        attributes(orientedMatrice.mtx)$interactions$orientation <- TRUE
        attributes(orientedMatrice.mtx)$interactions$submatrix.name <-
            attributes(orientedMatrice.mtx)$interactions$name
        names(orientedMatrice.mtx) <-
            attributes(orientedMatrice.mtx)$interactions$name
        return(orientedMatrice.mtx)
    } else {
        return(t(apply(
            as.data.frame(apply(matrice.mtx, 1, rev)),
            1,
            rev
        )))
    }
}
