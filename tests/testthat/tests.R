#==============================
# Libraries
#==============================
library(GenomicED)
library(GenomicRanges)
library(S4Vectors)

#==============================
# Data
#==============================
data("submatrixPF_Ctrl.mtx_lst")
data("anchors_Peaks.gnr")
data("baits_Peaks.gnr")
data("submatrixRF_Ctrl.mtx_lst")
data("submatrixRF.mtx_lst")
data("anchors_Index.gnr")
data("baits_Index.gnr")
data("interactions.gni")
data("HiC_ctrl.cmx_lst")
data("domains.gnr")

#==============================
# Global Variables
#==============================
seqlengths.num <- c('2L'=23513712, '2R'=25286936)
chromSize.dtf  <- data.frame(
  seqnames   = names(seqlengths.num ), 
  seqlengths = seqlengths.num
)
binSize.num <- 1000

#==============================
# Test ImportHiC
#==============================
# Download HiC File
options(timeout = 3600)
temp.dir <- file.path(tempdir(), "HIC_DATA")
dir.create(temp.dir)
Hic.url <- "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/7386f953-8da9-47b0-acb2-931cba810544/4DNFIOTPSS3L.hic"
HicOutput.pth <- file.path(temp.dir, "Control_HIC.hic")
download.file(Hic.url, HicOutput.pth, method = 'curl', extra = '-k')

Mcool.url <- "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/4f1479a2-4226-4163-ba99-837f2c8f4ac0/4DNFI8DRD739.mcool"
McoolOutput.pth <- file.path(temp.dir, "HeatShock_HIC.mcool")
download.file(Mcool.url, McoolOutput.pth, method = 'curl', extra = '-k')

# Import file in R
HiC_Ctrl.cmx_lst <- ImportHiC(
  file.pth    = HicOutput.pth,
  res.num     = binSize.num,
  chrom_1.chr = c("2L", "2R", "2L"),
  chrom_2.chr = c("2L", "2R", "2L")
)
HiC_HS.cmx_lst <- ImportHiC(
  file.pth    = McoolOutput.pth,
  res.num     = binSize.num,
  chrom_1.chr = c("2L", "2R", "2L"),
  chrom_2.chr = c("2L", "2R", "2L"),
  cores.num   = 2
  )
ImportHiC(
  file.pth    = HicOutput.pth,
  res.num     = "5Kb",
  chrom_1.chr = "2L"
)

#==============================
# Test NormalizeHiC
#==============================
HiC_Ctrl.cmx_lst <- NormalizeHiC(HiC_Ctrl.cmx_lst)
HiC_HS.cmx_lst <- NormalizeHiC(HiC_HS.cmx_lst)
NormalizeHiC(HiC_Ctrl.cmx_lst, interaction.type="cis", method.chr="VC", cores.num=2)
NormalizeHiC(HiC_Ctrl.cmx_lst, method.chr="trans", method.chr="VC_SQRT")
NormalizeHiC(HiC_Ctrl.cmx_lst, method.chr="all")

#==============================
# Test ExpectedHiC
#==============================
HiC_Ctrl.cmx_lst <- ExpectedHiC(HiC_Ctrl.cmx_lst)
HiC_HS.cmx_lst <-ExpectedHiC(HiC_HS.cmx_lst, cores.num=2)

#==============================
# Test IndexFeatures
#==============================
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

#==============================
# Test SearchPairs
#==============================
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

#==============================
# Test ExtractSubmatrix
#==============================
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

#==============================
# Test FilterInteractions
#==============================
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
  matrices.lst      = submatrixPF_Ctrl.mtx_lst,
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

#==============================
# Test GetQuantif
#==============================
GetQuantif(
  matrices.lst  = submatrixRF_Ctrl.mtx_lst,
  area.fun      = "center",
  operation.fun = "mean"
)
GetQuantif(
  matrices.lst  = submatrixPF_Ctrl.mtx_lst,
  area.fun      = "center",
  operation.fun = "mean"
)

#==============================
# Test OrienteMatrix
#==============================
OrienteMatrix(submatrixPF_Ctrl.mtx_lst)

#==============================
# Test Aggregation
#==============================
Aggregation(
  matrices.lst = submatrixPF_Ctrl.mtx_lst, 
  agg.fun      = "sum",
  trans.fun    = "qtl", 
  rm0.bln      = FALSE
)
diffAggreg.mtx <- Aggregation(
  ctrlMatrices.lst    = submatrixRF_Ctrl.mtx_lst,
  matrices.lst        = submatrixRF.mtx_lst,
  minDist             = "9Kb",
  maxDist             = "11Kb",
  agg.fun             = "mean",
  rm0.bln             = FALSE,
  diff.fun            = "substraction",
  scaleCorrection.bln = TRUE,
  correctionArea.lst  =  list(
    i = c(1:30),
    j = c(72:101)
    ),
  statCompare.bln = TRUE
)

#==============================
# Test ggAPA and PlotAPA
#==============================
PlotAPA(
    apa.mtx                  = diffAggreg.mtx,
    trimPrct.num             = 20,
    minBoundary.num          = -2,
    center.num               = 0,
    maxBoundary.num          = 2,
    minConditionBoundary.num = 0,
    maxConditionBoundary.num = 2
)
ggAPA(
    apa.mtx      = diffAggreg.mtx,
    title.chr    = "APA",
    center.num   = 0,
    trimPrct.num = NULL,
    bounds.chr   = "both",
    blurPass.num = 1,
    blurSd.num   = 0.5,
    heatmap.col  = viridis(6)
)

#==============================
# Complete Tests
#==============================
TrimOutliers(rnorm(1000))
GaussBox(scale.chr="int")
Gauss(x=1, y=1)
GenomicSystem(1000000000,2)
GRange_1.grn <- StrToGRanges(c("chr1:1-100:+","chr2:400-500:-"))
GRange_2.grn <- StrToGRanges("chr1:10-50:*")
MergeGRanges(list(GRange_1.grn,GRange_2.grn), reduce.bln=TRUE, sort.bln=TRUE)
Hue(paletteLength.num=1)
Hue(paletteLength.num=2)
PadMtx(mat.mtx=matrix(1:25,5,5),  padSize.num=1, value.num=NULL, side.chr=c('top','bot','right','left') )
ReduceRun(first.rle=rle(c("A","A","B")), second.rle=rle(c("A","B","B")), reduceFun.chr="paste", sep="_" )
