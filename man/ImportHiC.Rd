% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ImportHiC.R
\name{ImportHiC}
\alias{ImportHiC}
\title{Import Hic data}
\usage{
ImportHiC(
  file.pth = NULL,
  res.num = NULL,
  chromSize.dtf = NULL,
  chrom_1.chr = NULL,
  chrom_2.chr = NULL,
  verbose.bln = FALSE,
  cores.num = 1
)
}
\arguments{
\item{file.pth}{<GRanges or Pairs[GRanges] or GInteractions>: The genomic feature on which compute the extraction of HiC submatrix. Extension should be .hic, .cool, .mcool, .h5, .hdf5, .HDF5 or .bedpe" assuming .h5 et .hdf5 are only for cool (not mcool).}

\item{res.num}{<numeric>: The HiC resolution.}

\item{chromSize.dtf}{<data.frame>: A data.frame where first colum correspond to the chromosomes names, and the second column correspond to the chromosomes lengths in base pairs.}

\item{chrom_1.chr}{<numeric>: The seqnames of firsts chromosmes (rows in matrix).}

\item{chrom_2.chr}{<numeric>: The seqnames of second chromosmes (col in matrix). If is NULL is equal to chrom_1.chr (Defalt NULL)}

\item{verbose.bln}{<logical>: If TRUE show the progression in console. (Default FALSE)}

\item{cores.num}{<numerical> : An integer to specify the number of cores. (Default 1)}
}
\value{
A matrices list.
}
\description{
Import ..hic, .cool, .mcool or .bedpe data
}
\details{
ImportHiC
}
\examples{
\donttest{

# Prepare Temp Directory
options(timeout = 3600)
temp.dir <- file.path(tempdir(), "HIC_DATA")
dir.create(temp.dir)

# Download .hic file
Hic.url <- paste0(
    "https://4dn-open-data-public.s3.amazonaws.com/",
    "fourfront-webprod/wfoutput/",
    "7386f953-8da9-47b0-acb2-931cba810544/4DNFIOTPSS3L.hic"
)
HicOutput.pth <- file.path(temp.dir, "Control_HIC.hic")
download.file(Hic.url, HicOutput.pth, method = "curl", extra = "-k")

# Import .hic file
HiC_Ctrl.cmx_lst <- ImportHiC(
    file.pth = HicOutput.pth,
    res.num = 100000,
    chrom_1.chr = c("2L", "2L", "2R"),
    chrom_2.chr = c("2L", "2R", "2R")
)

# Download .mcool file
Mcool.url <- paste0(
    "https://4dn-open-data-public.s3.amazonaws.com/",
    "fourfront-webprod/wfoutput/",
    "4f1479a2-4226-4163-ba99-837f2c8f4ac0/4DNFI8DRD739.mcool"
)
McoolOutput.pth <- file.path(temp.dir, "HeatShock_HIC.mcool")
download.file(Mcool.url, McoolOutput.pth, method = "curl", extra = "-k")

# Import .mcool file
HiC_HS.cmx_lst <- ImportHiC(
    file.pth = McoolOutput.pth,
    res.num = 100000,
    chrom_1.chr = c("2L", "2L", "2R"),
    chrom_2.chr = c("2L", "2R", "2R")
)
}

}
