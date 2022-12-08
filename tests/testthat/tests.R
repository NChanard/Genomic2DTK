# Data
data("Beaf32_Peaks.gnr")
data("TSS_Peaks.gnr")
data("TADs_Domains.gnr")
data("HiC_Ctrl.cmx_lst")
data("HiC_HS.cmx_lst")

# Global Variables
seqlengths.num <- c('2L'=23513712, '2R'=25286936)
chromSize.dtf  <- data.frame(
    seqnames   = names(seqlengths.num ), 
    seqlengths = seqlengths.num
)

# Test BalanceHiC
BalanceHiC(HiC_Ctrl.cmx_lst, interaction.type=c("cis", "trans"), method.chr="ICE")
BalanceHiC(HiC_Ctrl.cmx_lst, interaction.type=c("cis", "trans"), method.chr="VC")
BalanceHiC(HiC_Ctrl.cmx_lst, interaction.type=c("cis", "trans"), method.chr="VC_SQRT")
BalanceHiC(HiC_Ctrl.cmx_lst, interaction.type="trans", method.chr="ICE")
BalanceHiC(HiC_Ctrl.cmx_lst, interaction.type="trans", method.chr="VC")
BalanceHiC(HiC_Ctrl.cmx_lst, interaction.type="trans", method.chr="VC_SQRT")
BalanceHiC(HiC_Ctrl.cmx_lst, interaction.type="cis", method.chr="ICE")
BalanceHiC(HiC_Ctrl.cmx_lst, interaction.type="cis", method.chr="VC")
BalanceHiC(HiC_Ctrl.cmx_lst, interaction.type="cis", method.chr="VC_SQRT")
BalanceHiC(HiC_Ctrl.cmx_lst, interaction.type="all", method.chr="ICE")
BalanceHiC(HiC_Ctrl.cmx_lst, interaction.type="all", method.chr="VC")
BalanceHiC(HiC_Ctrl.cmx_lst, interaction.type="all", method.chr="VC_SQRT")

HiC_Ctrl.cmx_lst <- BalanceHiC(HiC_Ctrl.cmx_lst, interaction.type="cis")
HiC_HS.cmx_lst   <- BalanceHiC(HiC_HS.cmx_lst,   interaction.type="cis")

# Test OverExpectedHiC
HiC_Ctrl.cmx_lst <- OverExpectedHiC(HiC_Ctrl.cmx_lst)
HiC_HS.cmx_lst   <- OverExpectedHiC(HiC_HS.cmx_lst, cores.num = 2)

# Test SwitchMatrix
SwitchMatrix(HiC_Ctrl.cmx_lst, matrixKind.chr="norm")
SwitchMatrix(HiC_Ctrl.cmx_lst, matrixKind.chr="o/e")

# Test IndexFeatures
Beaf32_Index.gnr <- IndexFeatures(
    gRange.gnr_lst        = list(Beaf=Beaf32_Peaks.gnr), 
    constraint.gnr        = TADs_Domains.gnr,
    chromSize.dtf         = chromSize.dtf,
    binSize.num           = 100000,
    method.chr            = "max",
    variablesName.chr_vec = "score",
    cores.num             = 2,
    verbose.bln           = TRUE
)
TSS_Index.gnr <- IndexFeatures(
    gRange.gnr_lst        = list(TSS=TSS_Peaks.gnr), 
    constraint.gnr        = TADs_Domains.gnr,
    chromSize.dtf         = chromSize.dtf,
    binSize.num           = 100000,
    method.chr            = "max",
    variablesName.chr_vec = "score",
    cores.num             = 2,
    verbose.bln           = TRUE
)
IndexFeatures(
    gRange.gnr_lst        = list(TSS_1=TSS_Peaks.gnr, TSS_2=TSS_Peaks.gnr), 
    constraint.gnr        = NULL,
    chromSize.dtf         = chromSize.dtf,
    binSize.num           = 100000,
    cores.num             = 1,
    verbose.bln           = TRUE
)
# Test SearchPairs
Beaf_TSS.gni <- SearchPairs(
    indexAnchor.gnr = Beaf32_Index.gnr,
    indexBait.gnr   = TSS_Index.gnr,
    minDist.num     = NULL, 
    maxDist.num     = NULL,
    cores.num       = 1,
    verbose.bln     = TRUE
)
SearchPairs(
    indexAnchor.gnr = Beaf32_Index.gnr,
    minDist.num     = "1", 
    maxDist.num     = "1MB",
    cores.num       = 1,
    verbose.bln     = TRUE
)
SearchPairs(
    indexAnchor.gnr = Beaf32_Index.gnr,
    indexBait.gnr   = TSS_Index.gnr,
    minDist.num     = 1, 
    maxDist.num     = 1000000,
    cores.num       = 2,
    verbose.bln     = TRUE
)


# Test ExtractSubmatrix
submatrixPF_Ctrl.mtx_lst <- ExtractSubmatrix(
    feature.gn         = Beaf_TSS.gni,
    hic.cmx_lst        = HiC_Ctrl.cmx_lst,
    res.num            = NULL,
    referencePoint.chr = "pf",
    matriceDim.num     = 21,
    cores.num          = 1,
    verbose.bln        = TRUE
)
submatrixPF_HS.mtx_lst <- ExtractSubmatrix(
    feature.gn         = Beaf_TSS.gni,
    hic.cmx_lst        = HiC_HS.cmx_lst,
    res.num            = NULL,
    referencePoint.chr = "pf",
    matriceDim.num     = 21,
    cores.num          = 1,
    verbose.bln        = TRUE
)
submatrixRF_Ctrl.mtx_lst <- ExtractSubmatrix(
    feature.gn         = Beaf_TSS.gni,
    hic.cmx_lst        = HiC_Ctrl.cmx_lst,
    res.num            = NULL,
    referencePoint.chr = "rf",
    matriceDim.num     = 21,
    cores.num          = 1,
    verbose.bln        = TRUE
)
ExtractSubmatrix(
    feature.gn         = TADs_Domains.gnr,
    hic.cmx_lst        = HiC_Ctrl.cmx_lst,
    referencePoint.chr = "rf",
    matriceDim.num     = 101,
    cores.num          = 2,
    verbose.bln        = FALSE
)


# Test FilterInteractions
target.lst <- list(
    anchor.Beaf.name = c("Beaf32_113"),
    bait.TSS.name    = c("FBgn0267378")
)
selection.fun = function(){
    Reduce(union, list(anchor.Beaf.name, bait.TSS.name) )
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

# Test GetQuantif
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

# Test OrientateMatrix
OrientateMatrix(submatrixPF_Ctrl.mtx_lst)

# Test Aggregation
Aggregation(
    matrices.lst = submatrixRF_Ctrl.mtx_lst, 
    agg.fun      = "sum",
    trans.fun    = "qtl", 
    rm0.bln      = FALSE
    )
diffAggreg.mtx <- Aggregation(
    ctrlMatrices.lst    = submatrixPF_Ctrl.mtx_lst,
    matrices.lst        = submatrixPF_HS.mtx_lst,
    minDist             = 1,
    maxDist             = "5Mb",
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

# Test ggAPA and PlotAPA
ggAPA(
    apa.mtx      = diffAggreg.mtx,
    title.chr    = "APA",
    colMin.num   = 0,
    colMax.num   = 10,
    trimPrct.num = 20,
    bounds.chr   = "both",
    blurPass.num = 1,
    blurSd.num   = 0.5,
    heatmap.col  = NULL
)


ggAPA(
    apa.mtx      = diffAggreg.mtx,
    title.chr    = "APA",
    colMid.num   = 0,
    trimPrct.num = NULL,
    bounds.chr   = "both",
    blurPass.num = 1,
    blurSd.num   = 0.5,
    heatmap.col  = viridis(6)
)

# Complete Tests
TrimOutliers(rnorm(1000))
GaussBox(scale.chr="int")
Gauss(x=1, y=1)
GenomicSystem(1000000000,2)
GRange_1.grn <- StrToGRanges(c("chr1:1-100:+","chr2:400-500:-"))
GRange_2.grn <- StrToGRanges("chr1:10-50:*")
MergeGRanges(list(GRange_1.grn,GRange_2.grn), reduce.bln=TRUE, sort.bln=TRUE)
Hue(paletteLength.num=1)
Hue(paletteLength.num=2)
PadMtx(mat.mtx=matrix(1:25,5,5), padSize.num=1, value.num=0, side.chr=c('top','bot','right','left') )
ReduceRun(first.rle=rle(c("A","A","B")), second.rle=rle(c("A","B","B")), reduceFun.chr="paste", sep="_" )
MeanScale(rnorm(500,500))
BreakVector(x.num=rnorm(500,500), n.num=50, method.chr="density")
BreakVector(x.num=rnorm(500,500), n.num=50, center.num=-500)
set.seed(123)
mat.spm = as(matrix(floor(runif(7*13,0,2)),7,13), "dgCMatrix")
Rise0(mat.spm=mat.spm, which.ndx=c(1,3,6,10,12))
Rise0(mat.spm=mat.spm, coord.dtf=data.frame(i=c(1,5,3), j=c(1,2,3) ) )
Rise0(mat.spm=mat.spm)
GetFileExtension(path.pth="my/path/to/my/file.txt")