suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))

## tile plot to show the log2-fold-change values of each SM cluster in its own TF
## overexpression dataset

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################

analysisName <- "self_cluster_DEG"
outDir <- here::here("analysis", "08_polII_analysis", "02_self_cluster_DEG")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

if(!dir.exists(outDir)){
  dir.create(outDir, recursive = T)
}

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
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF
col_geneId <- "GID"

col_lfc <- "log2FoldChange"
col_pval <- "pvalue"

###########################################################################

degIds <- suppressMessages(readr::read_tsv(file = file_degIds))

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% degIds$degId)


clusterInfo <- AnnotationDbi::select(
  x = orgDb, keys = rnaseqInfo$SM_TF, columns = c("GENE_NAME", "SM_CLUSTER", "SM_ID"), keytype = "GID"
) %>% 
  dplyr::filter(!is.na(SM_ID))


## choose only those SMTFs which are part of any SM cluster
rnaseqInfo <- dplyr::left_join(x = rnaseqInfo, y = clusterInfo, by = c("SM_TF" = "GID")) %>% 
  dplyr::filter(!is.na(SM_ID)) %>% 
  dplyr::mutate(
    geneLabel = paste(SM_TF, " (", GENE_NAME, ")", sep = ""),
    geneLabel = if_else(condition = SM_TF == GENE_NAME, true = SM_TF, false = geneLabel)
  )

smtfLables <- structure(rnaseqInfo$geneLabel, names = rnaseqInfo$SM_TF)

###########################################################################
## get DEG data for each SM_TF cluster from its own polII_DEG set
rowId <- 1
degData <- NULL


for (rowId in 1:nrow(rnaseqInfo)) {
  
  clusterGenes <- suppressMessages(
    AnnotationDbi::select(
      x = orgDb, keys = rnaseqInfo$SM_ID[rowId],
      columns = c("SM_CLUSTER", "GID", "GENE_NAME"), keytype = "SM_ID"
    )
  ) %>% 
    dplyr::mutate(
      SM_TF = rnaseqInfo$SM_TF[rowId],
      tfLabel = rnaseqInfo$geneLabel[rowId]
    ) %>% 
    dplyr::select(SM_TF, tfLabel, geneId = GID, GENE_NAME, everything())
  
  tmpDf <- suppressMessages(readr::read_tsv(file = rnaseqInfo$deg[rowId])) %>% 
    dplyr::select(geneId, !!col_lfc, shrinkLog2FC, pvalue, padj) %>% 
    dplyr::mutate(comparison = rnaseqInfo$comparison[rowId])
  
  clusterGenes <- dplyr::left_join(
    x = clusterGenes, y = tmpDf, by = "geneId"
  )
  
  degData <- dplyr::bind_rows(degData, clusterGenes)
  
}


# degData %>% 
#   dplyr::group_by(SM_TF) %>% 
#   dplyr::summarise(n = n()) %>% 
#   view()

genePos <- as.data.frame(genes(txDb)) %>% 
  dplyr::select(geneId = gene_id, chr = seqnames, start, end, strand)

degData <- dplyr::left_join(x = degData, y = genePos, by = "geneId") %>% 
  dplyr::arrange(SM_TF, chr, start)

readr::write_tsv(x = degData, file = paste(outPrefix, ".combined_DEGs.tab", sep = ""))

###########################################################################

plotData <- dplyr::group_by(degData, SM_TF) %>% 
  dplyr::arrange(chr, start, .by_group = TRUE) %>% 
  dplyr::mutate(index = row_number()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    significance = if_else(
      condition = !!sym(col_pval) <= cutoff_fdr, true = "significant",
      false = "non-significant", missing = "non-significant"
    )
  ) %>% 
  tidyr::unite(col = "tf_cluster_grp", tfLabel, SM_CLUSTER, sep = ": ")

pltTitle <- "RNA-polII ChIPseq log2(SMTF-OE / WT) of SM clusters in its own TF overexpression data"

pt_tiles <- ggplot() +
  geom_tile(
    data = plotData,
    mapping = aes(x = index, y = tf_cluster_grp, fill = shrinkLog2FC, color = significance),
    size = 0.2, height = 1) +
  scale_fill_gradient2(
    name = paste("log2(", "polII fold change", ")", sep = ""),
    low = "#313695", mid = "#F7F7F7", high = "#d73027", midpoint = 0,
    limit = c(-2.5, 2.5), oob = scales::squish
  ) +
  scale_colour_manual(
    name = "p-value <= 0.05",
    values = c("significant" = "black", "non-significant" = "white")
  ) +
  scale_x_continuous(expand = expansion(add = c(0.0, 0.0))) +
  facet_wrap(
    facets = . ~ tf_cluster_grp, scales = "free_y",
    ncol = 4, dir = "v"
  ) +
  ggtitle(pltTitle) +
  guides(color = guide_legend(override.aes = list(shape = 22, size = 2, fill = "grey"))) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.15, "lines"),
    strip.background = element_rect(fill="white", size = 0.2),
    strip.text.x = element_text(size = 12, hjust = 0),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    plot.margin = unit(rep(0.2, 4), "cm")
  )

# pt_tiles

ggsave(filename = paste(outPrefix, ".tiles.png", sep = ""), plot = pt_tiles, width = 14, height = 8)

###########################################################################

