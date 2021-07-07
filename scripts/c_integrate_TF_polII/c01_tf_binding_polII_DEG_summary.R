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
file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
file_peakStats <- here::here("analysis", "02_QC_TF", "TF_ChIP_summary.best_replicates.tab")
file_degStats <- here::here(
  "analysis", "08_polII_analysis", "01_polII_DEGs_summary", "polII_DEGs_summary.stats.tab"
)
file_SMTF_selfFoldChange <- here::here(
  "analysis", "08_polII_analysis", "01_polII_DEGs_summary", "SMTF_LFC_in_self_OE.tab"
)

TF_dataPath <- here::here("data", "TF_data")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

##################################################################################

productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::filter(has_polII_ChIP == "has_data", has_TF_ChIP == "has_data", copyNumber == "sCopy")

peakStats <- suppressMessages(readr::read_tsv(file = file_peakStats)) %>% 
  dplyr::select(sampleId, peaks_pval20)

degStats <- suppressMessages(readr::read_tsv(file = file_degStats)) %>% 
  dplyr::select(comparison, up, down, total)

selfLfc <- suppressMessages(readr::read_tsv(file = file_SMTF_selfFoldChange)) %>% 
  dplyr::select(degId, log2FoldChange, pvalue, padj, log10_padj, significant)

##################################################################################

dataSummary <- dplyr::left_join(x = productionData, y = peakStats, by = c("tfId" = "sampleId")) %>% 
  dplyr::left_join(y = degStats, by = c("degId" = "comparison")) %>% 
  dplyr::left_join(y = selfLfc, by = "degId") %>% 
  dplyr::mutate(degTotal = up + down) %>% 
  dplyr::arrange(peaks_pval20, degTotal) %>% 
  dplyr::mutate(
    SMTF = forcats::as_factor(SMTF),
    SMTF_name = forcats::as_factor(SMTF_name)
  )

readr::write_tsv(file = paste(outPrefix, ".stats.tab", sep = ""), x = dataSummary)

ptTheme <- theme_bw() +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 16),
    title = element_text(size = 18),
    panel.grid = element_blank()
  )

pt_peaks <- ggplot(data = dataSummary) +
  geom_bar(mapping = aes(y = SMTF_name, x = peaks_pval20),
           stat = "identity", fill = "black") +
  labs(title = "TF ChIPseq peak count") +
  ptTheme +
  theme(
    axis.title = element_blank()
  )


pt_deg <- dplyr::mutate(dataSummary, down = -1*down) %>% 
  ggplot(mapping = aes(y = SMTF_name)) +
  geom_bar(mapping = aes(x = down), stat = "identity", fill = "#313695") +
  geom_bar(mapping = aes(x = up), stat = "identity", fill = "#d73027") +
  labs(title = "OE/WT polII ChIPseq DEG count") +
  ptTheme

pt_merged <- ggarrange(pt_peaks, pt_deg, widths = c(1, 1), ncol = 2, align = "h")

png(filename = paste(outPrefix, ".stats.png", sep = ""), width = 4000, height = 3000, res = 300)
pt_merged
dev.off()

#############
## scatter plot of over-expression level (DEG-LFC) vs number of peak
pt_scatter <- ggplot(data = dataSummary, mapping = aes(x = peaks_pval20, y = log2FoldChange)) +
  geom_point(size = 4) +
  labs(
    title = "correlation between scale of over expression and number of peaks",
    y = "log2(OE/WT)", x = "# of peaks"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    title = element_text(size = 15),
    panel.grid = element_blank()
  )

ggsave(
  filename = paste(outPrefix, ".LFC_vs_peak_count_corr.pdf", sep = ""),
  plot = pt_scatter, width = 8, height = 8
)



