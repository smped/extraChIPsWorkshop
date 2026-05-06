#' This script was used to prepare the data for this workshop.
#' The full required data structure is:
#' data |
#'      |- annotations |
#'                     |- gencode.v44lift37.annotation.gtf.gz
#'                     |- hg19-blacklist.v2.bed.gz
#'                     |_ SRR8315192_greylist.bed
#'      |- bam  |
#'              |-
#'      |_ macs |
#'              |-
#'
#' These files contain a subset of data from the complete dataset
#' from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA509779&o=acc_s%3Aa
#' being processed at https://github.com/smped/PRJNA509779.
#' The original publication containing all relevant
#' biological background is at https://doi.org/10.1038/s41591-020-01168-7
#'
#' Two ChIP targets (ERa & H3K27ac) were used in ZR-75-1 cells under control
#' conditions (E2), followed by the addition of dihydro testosterone (E2+DHT).
#'
#' Currently all files contain data for the complete genome but may need to be
#' stripped back to a single chromosome for a lighter workshop.
#' If possible, try using all data saved as counts instead of using the bam
#' files.

## Use cached data to ensure it's only downloaded once
library(BiocFileCache)
bfc <- BiocFileCache(file.path("data"), ask = FALSE)

## Use a helper function to download each folder & unzip
download_and_unzip <- function(bfc, url, destdir = "data") {
  ## Only downloads if not already cached
  path <- bfcrpath(bfc, url)
  message("Extracting ", basename(url), " ...")
  unzip(path, exdir = destdir)
  invisible(destdir)
}

## The annotation files
download_and_unzip(
  bfc, "https://zenodo.org/records/10052122/files/annotations.zip"
)

## The bam files. NB: these may take 30mins or so
download_and_unzip(bfc, "https://zenodo.org/records/10052122/files/bam.zip")

## The macs2 output
download_and_unzip(bfc, "https://zenodo.org/records/10052122/files/macs2.zip")
