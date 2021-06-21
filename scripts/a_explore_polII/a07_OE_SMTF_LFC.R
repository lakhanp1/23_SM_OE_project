suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))

## Figure to show the log2(OE/WT) for SMTF in its own overexpression polII-ChIPseq data

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################

analysisName <- "SM_OE_diff_summary"
outDir <- here::here("analysis", "08_polII_analysis", "01_polII_DEGs_summary")
outPrefix <- paste(outDir, "/", sep = "")


file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
diffDataPath <- here::here("analysis", "06_polII_diff")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_degStats <- here::here(
  "analysis", "08_polII_analysis", "01_polII_DEGs_summary", "polII_DEGs_summary.stats.tab"
)

useAllGroupsSamples <- FALSE

cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1

orgDb <- org.Anidulans.FGSCA4.eg.db
col_geneId <- "GID"

###########################################################################
productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::filter(has_polII_ChIP == "has_data", has_TF_ChIP == "has_data", copyNumber == "sCopy")

degStats <- suppressMessages(readr::read_tsv(file = file_degStats))

productionData <- dplyr::left_join(
  x = productionData, y = degStats, by = c("degId" = "comparison")
)

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% productionData$degId)

rnaseqInfoList <- purrr::transpose(rnaseqInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "comparison"))

productionData$log2FoldChange <- NA
productionData$pvalue <- NA
productionData$padj <- NA


rowId <- 1

for (rowId in 1:nrow(productionData)) {
  
  degInfo <- purrr::pluck(.x = rnaseqInfoList, productionData$degId[[rowId]])
  
  diffData <- suppressMessages(readr::read_tsv(file = degInfo$deseq2)) %>% 
    dplyr::filter(geneId == degInfo$SM_TF)
  
  productionData$log2FoldChange[rowId] <- diffData$log2FoldChange
  productionData$pvalue[rowId] <- diffData$pvalue
  productionData$padj[rowId] <- diffData$padj
  
}


plotData <-  productionData %>% 
  dplyr::arrange(log2FoldChange) %>% 
  dplyr::mutate(
    log10_padj = -log10(padj),
    significant = if_else(condition = padj <= 0.05, "significant", "non-significant"),
    SM_TF = forcats::as_factor(SM_TF),
    geneName = forcats::as_factor(geneName)
  )


###########################################################################

## color scales
# scaleLim <- ceiling(min(5, max(goData[[logPvalCol]])))
scaleLim <- 5
brk = c(1, 1.30103, 2:scaleLim, scaleLim+1)
scaleLabels <- c(format(1/(10^c(1, 1.30103, 2:scaleLim)), drop0trailing = T, scientific = T), "smaller")



pt_lfc <- ggplot(
  data = plotData,
  mapping = aes(x = log2FoldChange, y = geneName)) +
  geom_point(mapping = aes(fill = log10_padj, color = significant),
             shape = 21, size = 5, stroke = 1) +
  scale_fill_gradientn(
    name = "-log10(padj)",
    values = scales::rescale(c(seq(1, max(plotData$log10_padj), length.out = 9))),
    colours = RColorBrewer::brewer.pal(n = 9, name = "PuRd"),
    breaks = brk,
    labels = scaleLabels,
    guide = guide_colorbar(barheight = 10, draw.llim = FALSE, order = 1),
    oob = squish,
    limits = c(1, scaleLim + 1)
  ) +
  scale_color_manual(
    name = "p-adjusted <= 0.05",
    values = c("significant" = "black", "non-significant" = "red")
  ) +
  labs(
    title = stringr::str_wrap("log2(fold-change) for SM TF genes in its own overexpression polII ChIPseq data"),
    x = "log2(fold-change)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    title = element_text(size =14, face = "bold"),
    legend.text = element_text(size = 12)
  )

ggsave(filename = paste(outPrefix, "SMTF_LFC_in_self_OE.pdf", sep = ""), plot = pt_lfc, width = 10, height = 10)
ggsave(filename = paste(outPrefix, "SMTF_LFC_in_self_OE.png", sep = ""), plot = pt_lfc, width = 10, height = 8, scale = 1.1)


## scatter plot to show correlation between SMTF overexpression level and number of DEGs
pt_scatter <- ggplot(data = plotData, mapping = aes(x = total, y = log2FoldChange)) +
  geom_point(size = 3) +
  geom_smooth(method=lm, se = FALSE, formula = y ~ x, color = "red") +
  ggpubr::stat_cor(method = "pearson", size = 6, label.x.npc = 0.5) +
  labs(
    x = "# of DEGs",
    y = "log2(OE / WT) for SMTFs",
    title = "Correlation between the level of over-expression and number of DEGs"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0),
    panel.grid = element_blank(),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    title = element_text(size =14, face = "bold"),
    legend.text = element_text(size = 12)
  )

ggsave(
  filename = paste(outPrefix, "SMTF_LFC_vs_DEGs_count.pdf", sep = ""),
  plot = pt_scatter, width = 9, height = 9
)





