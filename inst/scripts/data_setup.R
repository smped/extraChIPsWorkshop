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
library(tidyverse)
library(Rsamtools)
library(rtracklayer)
library(extraChIPs)
library(plyranges)
data_path <- file.path("~", "extraChIPs_data")
if (!dir.exists(data_path)) dir.create(data_path)
bfc <- BiocFileCache(data_path, ask = FALSE)

#### Download all/any that are required ####
url_list <- file.path(
  "https://zenodo.org/records/10052122/files",
  c("annotations.zip", "bam.zip", "macs2.zip")
)
lapply(url_list, \(x) bfcadd(bfc, x))

## Extract as zip files
lapply(bfcinfo(bfc)$rpath, unzip, exdir = data_path)

#### Index Bam Files ####
## The BamFiles may have conflicting timestamps with the indexes
## Reindex just for convenience
file.path(data_path, "bam") |>
  list.files(pattern = "bam$", full.names = TRUE) |>
  lapply(indexBam)

#### Setup the annotation files ####
extdata <- here::here("inst", "extdata", "annotation")
if (!dir.exists(extdata)) dir.create(extdata, recursive = TRUE)
## Copy the black & grey lists
file.copy(
  file.path(data_path, "annotations", "hg19-blacklist.v2.bed.gz"),
  file.path(extdata, "hg19-blacklist.v2.bed.gz"), overwrite = TRUE
)
file.copy(
  file.path(data_path, "annotations", "SRR8315192_greylist.bed"),
  file.path(extdata, "SRR8315192_greylist.bed"), overwrite = TRUE
)
sq <- defineSeqinfo("GRCh37", chr = TRUE)
gtf <- file.path(data_path, "annotations", "gencode.v44lift37.annotation.gtf.gz") |>
  import.gff(feature.type = c("gene", "transcript", "exon")) |>
  subset(seqnames %in% seqlevels(sq))
f <- gzfile(file.path(extdata, "gencode.v44lift37.reduced.gtf.gz"))
gtf |>
  select(
    source, type, score, phase,
    starts_with("gene"), starts_with("transcript"), starts_with("exon"),
    -contains(c("status", "support"))
  ) |>
  write_gff(f)
close(f)

#### Copy the peaks across ####
extdata <- here::here("inst", "extdata", "peaks")
if (!dir.exists(extdata)) dir.create(extdata, recursive = TRUE)
file.path(data_path, "macs2") |>
  list.files(full.names = TRUE, pattern = "narrowPeak$", recursive = TRUE) |>
  lapply(
    \(x) {
      ## Compress for convenience
      dest <- file.path(extdata, basename(x)) |>
        paste0(".gz") |>
        str_replace("E2", "H3K27ac_E2") |>
        gzfile(open = "w")
      x |>
        readLines() |>
        write(dest)
      close(dest)
    }
  ) |>
  invisible()

#### Prepare the BigWig objects ####
gr <- GRanges("chr22:30760000-30830000")
bwfl <- list.files(
  data_path, pattern = ".bw$", recursive = TRUE, full.names = TRUE
) %>%
  BigWigFileList()
names(bwfl) <- bwfl %>%
  path() %>%
  str_remove_all(".+macs2/") %>%
  str_replace_all("/", "_") %>%
  str_remove("_treat")
cov <- lapply(bwfl, import.bw, which = gr)
## And the export itself
extdata <- here::here("inst", "extdata", "bigwig")
if (!dir.exists(extdata)) dir.create(extdata, recursive = TRUE)
for (i in seq_along(cov)) {
  f <- file.path(extdata, names(cov)[[i]])
  export.bw(cov[[i]], f)
}

