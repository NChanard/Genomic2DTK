library(GenomicED)
library(GenomicRanges)
library(S4Vectors)

data("anchors_Peaks.gnr")
data("baits_Peaks.gnr")
seqlengths.num <- c('2L'=23513712, '2R'=25286936)
chromSize.dtf  <- data.frame(
  seqnames   = names(seqlengths.num ), 
  seqlengths = seqlengths.num
  )
binSize.num <- 1000
IndexFeatures(
  gRange.gnr_lst        = list(Beaf=anchors_Peaks.gnr, TSS=baits_Peaks.gnr), 
  constraint.gnr        = domains.gnr,
  chromSize.dtf         = chromSize.dtf,
  binSize.num           = binSize.num,
  method.chr            = "max",
  variablesName.chr_vec = "score",
  cores.num             = 1,
  verbose.bln           = TRUE
  )
IndexFeatures(
  gRange.gnr_lst        = list(Beaf=anchors_Peaks.gnr, TSS=baits_Peaks.gnr), 
  constraint.gnr        = NULL,
  chromSize.dtf         = chromSize.dtf,
  binSize.num           = binSize.num,
  method.chr            = "max",
  variablesName.chr_vec = "score",
  cores.num             = 2,
  verbose.bln           = FALSE
  )


data("anchors_Index.gnr")
data("baits_Index.gnr")
SearchPairs(
    indexAnchor.gnr = anchors_Index.gnr,
    indexBait.gnr   = baits_Index.gnr,
    minDist.num     = NULL, 
    maxDist.num     = NULL,
    cores.num       = 2,
    verbose.bln     = FALSE
    )
SearchPairs(
    indexAnchor.gnr = anchors_Index.gnr,
    minDist.num     = "1", 
    maxDist.num     = "50Kb",
    cores.num       = 1,
    verbose.bln     = TRUE
    )
SearchPairs(
    indexAnchor.gnr = anchors_Index.gnr,
    indexBait.gnr   = baits_Index.gnr,
    minDist.num     = 1, 
    maxDist.num     = 50000,
    cores.num       = 2,
    verbose.bln     = FALSE
    )
data("interactions.gni")
data("HiC_ctrl.cmx_lst")
data("domains.gnr")
ExtractSubmatrix(
  feature.gn         = interactions.gni,
  hic.cmx_lst        = HiC_ctrl.cmx_lst,
  res.num            = NULL,
  referencePoint.chr = "rf",
  matriceDim.num     = 101,
  cores.num          = 1,
  verbose.bln        = TRUE
  )

ExtractSubmatrix(
  feature.gn         = interactions.gni,
  hic.cmx_lst        = HiC_ctrl.cmx_lst,
  res.num            = NULL,
  referencePoint.chr = "pf",
  matriceDim.num     = 101,
  cores.num          = 2,
  verbose.bln        = TRUE
  )
ExtractSubmatrix(
  feature.gn         = domains.gnr,
  hic.cmx_lst        = HiC_ctrl.cmx_lst,
  referencePoint.chr = "rf",
  matriceDim.num     = 101,
  cores.num          = 2,
  verbose.bln        = FALSE
  )

data("submatrixPF_Ctrl.mtx_lst")
target.lst <- list(
  anchor.Beaf.name = c("Beaf32_8","Beaf32_15"),
  bait.Tss.name    = c("FBgn0031214","FBgn0005278"),
  name             = c("2L:74_2L:77"),
  distance         = function(columnElement){
    return(14000==columnElement || columnElement == 3000)
    }
  )

selection.fun = function(){
  Reduce(intersect, list(anchor.Beaf.name, bait.Tss.name ,distance) ) |>
  setdiff(name)
  }

FilterInteractions(
  matrices.lst = submatrixPF_Ctrl.mtx_lst,
  target.lst        = target.lst,
  selection.fun     = selection.fun
  )

FilterInteractions(
  interarctions.gni = attributes(submatrixPF_Ctrl.mtx_lst)$interactions,
  target.lst        = target.lst,
  selection.fun     = NULL
  )

target.lst <- list(interactions = attributes(submatrixPF_Ctrl.mtx_lst)$interactions[1:2])
FilterInteractions(
  interarctions.gni = attributes(submatrixPF_Ctrl.mtx_lst)$interactions,
  target.lst        = target.lst,
  selection.fun     = NULL
  )
 
target.lst <- list(first = InteractionSet::anchors(attributes(submatrixPF_Ctrl.mtx_lst)$interactions)[["first"]][1:2])
FilterInteractions(
  interarctions.gni = attributes(submatrixPF_Ctrl.mtx_lst)$interactions,
  target.lst        = target.lst,
  selection.fun     = NULL
  )

OrienteMatrix(submatrixPF_Ctrl.mtx_lst)
GetQuantif(
  matrices.lst  = submatrixPF_Ctrl.mtx_lst,
  area.fun      = "center",
  operation.fun = "mean"
  )

data("submatrixRF_Ctrl.mtx_lst")
GetQuantif(
  matrices.lst  = submatrixRF_Ctrl.mtx_lst,
  area.fun      = "center",
  operation.fun = "mean"
  )
Aggregation(
  matrices.lst = submatrixPF_Ctrl.mtx_lst, 
  agg.fun      = "sum",
  rm0.bln      = FALSE
  )
Aggregation(
  matrices.lst = submatrixPF_Ctrl.mtx_lst, 
  agg.fun      = "sum",
  rm0.bln      = TRUE,
  minDist      = 9000,
  maxDist      = 11000
  )

diffAggreg.mtx <- Aggregation(
  ctrlMatrices.lst    = submatrixRF_Ctrl.mtx_lst,
  matrices.lst        = submatrixRF.mtx_lst,
  minDist             = NULL,
  maxDist             = NULL,
  agg.fun             = "mean",
  rm0.bln             = FALSE,
  diff.fun            = "substraction",
  scaleCorrection.bln = TRUE,
  correctionArea.lst  =  list(
    i = c(1:30),
    j = c(72:101)
    ),
  statCompare.bln = TRUE)

PlotAPA(
    apa.mtx                  = diffAggreg.mtx,
    trimPrct.num             = 20,
    minBoundary.num          = -2,
    center.num               = 0,
    maxBoundary.num          = 2,
    minConditionBoundary.num = 0,
    maxConditionBoundary.num = 2
)