#' Work through the complete workflow here. It's based on the half-written
#' https://github.com/smped/extraChIPs_paper/blob/main/manuscript.Rmd, which
#' was abandoned when F1000R changed their rules & broke the Bioconductor
#' pipeline.
#'
#'


library(tidyverse)
library(extraChIPs)
library(rtracklayer)
library(plyranges)
library(scales)
library(glue)
library(ComplexUpset)
library(Rsamtools)
library(csaw)
library(BiocParallel)
library(quantro)
library(ggrepel)
library(patchwork)
library(xtable)
library(ggside)
register(MulticoreParam(4))
theme_set(theme_bw())

##### ERa: Fixed Width Peaks ######
## Define the seqinfo
sq <- defineSeqinfo("GRCh37", chr = TRUE)

## Define the samples data frame
samples_df <- tibble(
  id = paste0("SRR831518", seq(0, 5)),
  target = "ERa",
  treat = rep(c("E2", "E2+DHT"), each = 3) |> fct_inorder(),
  file = paste0(id, "_peaks.narrowPeak")
)
peak_path <- here::here("data", "macs2", samples_df$id, samples_df$file)

## Start by loading up the ERa peaks then forming the consensus peaks
peaks <- importPeaks(peak_path, centre = TRUE, seqinfo = sq)
names(peaks) <- names(peaks) %>% str_remove_all("_peaks.+")
peaks

treat_colours <- setNames(
  RColorBrewer::brewer.pal(9, "Set1")[2:1], levels(samples_df$treat)
)
peaks %>%
  map_int(length) %>%
  enframe(name = "id", value = "n_peaks") %>%
  left_join(samples_df) %>%
  ggplot(aes(id, n_peaks, fill = treat)) +
  geom_col() +
  scale_fill_manual(values = treat_colours) +
  scale_y_continuous(expand = expansion(c(0, 0.05)), labels = comma)

plotGrlCol(
  peaks, var = "width", geom = "violin", q = 0.01,
  df = samples_df, .id = "id", fill = "treat",
  trim = FALSE, draw_quantiles = 0.5, width = 0.8
) +
  labs(x = "ID", y = "Peak Width (bp)", fill = "Treatment") +
  scale_x_discrete(labels = \(x) str_remove(x, "_peaks.+")) +
  scale_y_log10() +
  scale_fill_manual(values = treat_colours)

plotGrlCol(
  peaks, var = "score", geom = "boxplot", q = 0.5, total_geom = "none",
  df = samples_df, .id = "id", fill = "treat"
) +
  labs(x = "ID", y = "Peak Width (bp)", fill = "Treatment") +
  scale_x_discrete(labels = \(x) str_remove(x, "_peaks.+")) +
  scale_y_log10() +
  scale_fill_manual(values = treat_colours)

plotOverlaps(peaks, sort_sets = "none", set_col = rep(treat_colours, each = 3))

## Consensus Peaks
consensus_by_treat <- peaks %>%
  split(samples_df$treat) %>%
  lapply(
    makeConsensus,
    p = 2/3, var = "centre", method = "coverage", min_width = 150
  ) %>%
  lapply(
    mutate,
    centre = map_int(centre, \(x) as.integer(mean(x)))
  )

consensus_by_treat %>%
  GRangesList() %>%
  plotOverlaps(set_col = treat_colours)


union_peaks <- consensus_by_treat %>%
  GRangesList() %>%
  makeConsensus(var = "centre") %>%
  mutate(centre = floor(map_dbl(centre, mean)))
