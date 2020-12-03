suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggpubr))


## generate TF summary plot heatmap
## TF ChIPseq correlation heatmap +
## TF peak count bar chart +
## peak occupancy heat map

rm(list = ls())

##################################################################################

analysisName <- "TF_polII_integration"
outDir <- here::here("analysis", "10_TF_polII_integration")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_dataSummary <- here::here("data", "reference_data", "raw_data_summary.tab")
file_peakSummary <- here::here("analysis", "02_QC_TF", "TF_ChIP_summary.best_replicates.tab")
file_degSummary <- here::here("analysis", "06_polII_diff", "polII_DEG.stats.tab")

TF_dataPath <- here::here("data", "TF_data")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

##################################################################################

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

dataSummary <- suppressMessages(readr::read_tsv(file = file_dataSummary)) %>% 
  dplyr::filter(has_TF_ChIP == "has_data", has_polII_ChIP == "has_data")

peakSummary <- suppressMessages(readr::read_tsv(file = file_peakSummary)) %>% 
  dplyr::select(sampleId, peaks_pval20)

degSummary <- suppressMessages(readr::read_tsv(file = file_degSummary))

##################################################################################

dataSummary <- dplyr::left_join(x = dataSummary, y = peakSummary, by = c("tfId" = "sampleId")) %>% 
  dplyr::left_join(y = degSummary, by = c("degId" = "comparison")) %>% 
  dplyr::mutate(degTotal = up + down) %>% 
  dplyr::arrange(peaks_pval20, degTotal) %>% 
  dplyr::mutate(
    geneId = forcats::as_factor(geneId)
  )

ptTheme <- theme_bw() +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 16),
    axis.title = element_blank(),
    title = element_text(size = 18),
    panel.grid = element_blank()
  )

pt_peaks <- ggplot(data = dataSummary) +
  geom_bar(mapping = aes(y = geneId, x = peaks_pval20),
           stat = "identity", fill = "black") +
  labs(title = "TFOE ChIPseq peak count") +
  ptTheme


pt_deg <- dplyr::mutate(dataSummary, down = -1*down) %>% 
  ggplot(mapping = aes(y = geneId)) +
  geom_bar(mapping = aes(x = down), stat = "identity", fill = "#313695") +
  geom_bar(mapping = aes(x = up), stat = "identity", fill = "#d73027") +
  labs(title = "TFOE/WT polII ChIPseq DEG count") +
  ptTheme

pt_merged <- ggarrange(pt_peaks, pt_deg, widths = c(1, 1), ncol = 2, align = "h")

png(filename = paste(outPrefix, ".stats.png", sep = ""), width = 4000, height = 3000, res = 300)
pt_merged
dev.off()








