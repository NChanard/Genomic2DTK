NotIn("A", c("A","B","C"))
NotIn("A", c("B","C","D"))
 
start.tim <- Sys.time()
for(i in seq_len(10000)){
 ShowLoading(start.tim, i , 10000)
}

start.tim <- Sys.time()
for(i in seq_len(10000)){
    Sys.sleep(3)
    ShowLoading(start.tim, i , 10000)
    if (i==3){
        break
    }
}

df1 <- data.frame(a = c(1:5), b = c(6:10))
df2 <- data.frame(a = c(11:15), b = c(16:20), c = LETTERS[1:5])
BindFillRows(df1,df2)
BindFillRows(list(df1,df2))

GRange.gnr <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(c("chr1", "chr2"), c(3, 1)),
    ranges = IRanges::IRanges(c(1,201,251,1), end = c(200,250,330,100), names = letters[1:4]),
    strand = S4Vectors::Rle(BiocGenerics::strand(c("*")), 4),
    score = c(50,NA,100,30)
    )
GRange.gnr
chromSize.dtf = data.frame(c("chr1","chr2"),c(350,100))
binSize.num <- 100
binnedGRanges.gnr <- BinGRanges(
    gRange.gnr = GRange.gnr,
    chromSize.dtf=chromSize.dtf,
    binSize.num=binSize.num,
    method.chr ="mean",
    variablesName.chr_vec="score",
    na.rm=TRUE
)
binnedGRanges.gnr <- BinGRanges(
    gRange.gnr = GRange.gnr,
    chromSize.dtf=chromSize.dtf,
    binSize.num=binSize.num,
    method.chr ="mean",
    variablesName.chr_vec="score",
    na.rm=TRUE,
    cores.num=2
)

StrToGRanges("chr1:1-100:+")
StrToGRanges(c("chr1:1-100:+","chr2:400-500:-","chr1:10-50:*"))

GenomicSystem(1540,3)
GenomicSystem(1540,2)
GenomicSystem(10,2)
GenomicSystem(1000,2)
GenomicSystem(1000000,2)
GenomicSystem(1000000000,2)
GenomicSystem("1Gbp")
GenomicSystem("1Mbp")
GenomicSystem("1Kbp")
GenomicSystem("10Bp")

GRange_1.grn <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(c("chr1", "chr2", "chr1"), c(1, 3, 1)),
    ranges = IRanges::IRanges(101:105, end = 111:115, names = letters[1:5]),
    strand = S4Vectors::Rle(BiocGenerics::strand(c("-", "+", "*", "+")), c(1, 1, 2, 1)),
    score = 1:5
)
GRange_2.grn <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(c("chr1", "chr3"), c(1, 4)),
    ranges = IRanges::IRanges(106:110, end = 116:120, names = letters[6:10]),
    strand = S4Vectors::Rle(BiocGenerics::strand(c("*", "+", "-")), c(2, 1, 2)),
    score = 6:10
)
MergeGRanges(GRange_1.grn,GRange_2.grn)
GRange.lst = list(GRange_1.grn,GRange_2.grn)
MergeGRanges(GRange.lst)
MergeGRanges(GRange.lst, reduce.bln=TRUE, sort.bln=TRUE)

GRange.grn <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(c("chr1", "chr2", "chr1"), c(1, 3, 1)),
    ranges = IRanges::IRanges(101:105, end = 111:115, names = letters[1:5]),
    strand = S4Vectors::Rle(BiocGenerics::strand(c("-", "+", "*", "+")), c(1, 1, 2, 1)),
    seqinfo = c(chr1=200, chr2=300),
    score = 1:5
)
SeqEnds(GRange.grn)

StrToGRanges("chr1:1-100:+")
StrToGRanges(c("chr1:1-100:+","chr2:400-500:-","chr1:10-50:*"))



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
  trans.fun    = "qtl", 
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