suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))



rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################

analysisName <- "respiration_DEG"
outDir <- here::here("analysis", "08_polII_diff_downstream", "05_mitochondria")
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
  x = orgDb, keys = rnaseqInfo$SM_TF, columns = c("GENE_NAME", "SM_CLUSTER"), keytype = "GID"
) %>% 
  dplyr::group_by(GID) %>% 
  dplyr::summarise(
    geneName = unique(GENE_NAME),
    SM_cluster = paste(unique(SM_CLUSTER), collapse = "/")
  )

rnaseqInfo <- dplyr::left_join(x = rnaseqInfo, y = clusterInfo, by = c("SM_TF" = "GID")) %>% 
  dplyr::mutate(
    geneLabel = paste(SM_TF, " (", geneName, ")", sep = ""),
    geneLabel = if_else(condition = SM_TF == geneName, true = SM_TF, false = geneLabel),
    degLabel = paste(geneLabel, ": ", SM_cluster, sep = "")
  )

## extract genes belonging to "GO:0022900 electron transport chain"
geneset <- AnnotationDbi::select(
  x = orgDb, keys = "GO:0022900", columns = c("GID", "GENE_NAME"), keytype = "GOALL"
) %>% 
  dplyr::select(geneId = GID, GENE_NAME) %>% 
  dplyr::mutate(
    geneLabel = paste(geneId, "(", GENE_NAME, ")", sep = ""),
    geneLabel = if_else(condition = geneId == GENE_NAME, true = geneId, false = geneLabel)
  )

genePos <- as.data.frame(genes(txDb)) %>% 
  dplyr::select(geneId = gene_id, chr = seqnames, start, end, strand)

geneset <- dplyr::left_join(x = geneset, y = genePos, by = "geneId") %>% 
  dplyr::arrange(chr, start) %>% 
  dplyr::mutate(
    geneId = forcats::as_factor(geneId)
  )

###########################################################################
## get DEG data for each SM_TF cluster from its own polII_DEG set
rowId <- 1
degData <- NULL


for (rowId in 1:nrow(rnaseqInfo)) {
  
  tmpDf <- suppressMessages(readr::read_tsv(file = rnaseqInfo$deg[rowId])) %>% 
    dplyr::select(geneId, !!col_lfc, shrinkLog2FC, pvalue, padj) %>% 
    dplyr::mutate(
      comparison = rnaseqInfo$comparison[rowId],
      SM_TF = rnaseqInfo$SM_TF[rowId],
      degLabel = rnaseqInfo$degLabel[rowId]
      )
  
  subData <- dplyr::left_join(
    x = geneset, y = tmpDf, by = "geneId"
  )
  
  degData <- dplyr::bind_rows(degData, subData)
  
}

readr::write_tsv(x = degData, file = paste(outPrefix, ".combined_DEGs.tab", sep = ""))

###########################################################################


lfcMat <- tidyr::pivot_wider(
  data = degData,
  id_cols = c(SM_TF),
  names_from = geneId,
  values_from = shrinkLog2FC
) %>% 
  tibble::column_to_rownames(var = "SM_TF") %>% 
  dplyr::select(!!!geneset$geneId) %>% 
  as.matrix()


fcHeatmap <- Heatmap(
  matrix = lfcMat,
  col = colorRamp2(breaks = c(-3, 0, 3), colors = c("#313695", "gray", "#d73027"), space = "LAB"),
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  cluster_columns = FALSE,
  show_row_names = TRUE,
  column_labels = geneset$geneLabel,
  row_labels = rnaseqInfo$degLabel,
  row_names_gp = gpar(fontsize = 16),
  row_names_max_width = unit(12, "cm"),
  column_names_gp = gpar(fontsize = 16), 
  # width = unit(4, "cm"),
  heatmap_legend_param = list(title = "\nlog2(OE/WT)")
)


pdf(file = paste(outPrefix, ".pdf", sep = ""), height = 11, width = 15)

draw(
  object = fcHeatmap,
  column_title = "RNA polII ChIPseq log2(OE/WT) for genes related to electron transport chain",
  column_title_gp = gpar(fontsize = 20, fontface = "bold")
)

dev.off()



