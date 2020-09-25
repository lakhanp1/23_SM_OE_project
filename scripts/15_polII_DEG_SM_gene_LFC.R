suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))

## plot fold change of SM gene in its own OE strain

rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################

analysisName <- "SM_OE_diff_summary"
outDir <- here::here("analysis", "08_polII_diff_downstream")
outPrefix <- paste(outDir, "/", analysisName, sep = "")


diffDataPath <- here::here("analysis", "06_polII_diff")
file_sampleInfo <- here::here("data", "reference_data", "polII_sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_degIds <- here::here("data", "reference_data", "production_data.polII_DEG_ids.txt")

useAllGroupsSamples <- FALSE

cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1

orgDb <- org.Anidulans.FGSCA4.eg.db
col_geneId <- "GID"


###########################################################################

degIds <- suppressMessages(readr::read_tsv(file = file_degIds))

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% degIds$comparison)

rnaseqInfo$log2FoldChange <- NA
rnaseqInfo$pvalue <- NA
rnaseqInfo$padj <- NA


row <- 1
for (row in 1:nrow(rnaseqInfo)) {
  diffData <- suppressMessages(readr::read_tsv(file = rnaseqInfo$deseq2[row])) %>% 
    dplyr::filter(geneId == rnaseqInfo$SM_TF[row])
  
  rnaseqInfo$log2FoldChange[row] <- diffData$log2FoldChange
  rnaseqInfo$pvalue[row] <- diffData$pvalue
  rnaseqInfo$padj[row] <- diffData$padj
  
}

plotData <-  rnaseqInfo %>% 
  dplyr::arrange(log2FoldChange) %>% 
  dplyr::mutate(
    log10_padj = -log10(padj),
    significant = if_else(condition = padj <= 0.05, "significant", "non-significant"),
    SM_TF = forcats::as_factor(SM_TF)
  )



## color scales
# scaleLim <- ceiling(min(5, max(goData[[logPvalCol]])))
scaleLim <- 5
brk = c(1, 1.30103, 2:scaleLim, scaleLim+1)
scaleLabels <- c(format(1/(10^c(1, 1.30103, 2:scaleLim)), drop0trailing = T, scientific = T), "smaller")



pt <- ggplot(
  data = plotData,
  mapping = aes(x = log2FoldChange, y = SM_TF)) +
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

ggsave(filename = paste(outPrefix, ".SM_gene_LFC.pdf", sep = ""), plot = pt, width = 10, height = 10)
ggsave(filename = paste(outPrefix, ".SM_gene_LFC.png", sep = ""), plot = pt, width = 10, height = 10, scale = 1.1)









