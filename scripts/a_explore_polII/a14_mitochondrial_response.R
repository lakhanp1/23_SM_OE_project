suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))



rm(list = ls())

source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################

analysisName <- "respiration_DEG"
outDir <- here::here("analysis", "08_polII_analysis", "mitochondria")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
diffDataPath <- here::here("analysis", "06_polII_diff")
file_sampleInfo <- here::here("data", "reference_data", "polII_sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_polIIDegs <- here::here("analysis", "06_polII_diff", "polII_DEGs.fpkm_filtered.combined.tab")

cutoff_fpkm <- 10
cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

col_geneId <- "GID"
col_lfc <- "log2FoldChange"
col_pval <- "padj"

###########################################################################

if(!dir.exists(outDir)){
  dir.create(outDir, recursive = T)
}

productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::filter(has_polII_ChIP == "has_data", has_TF_ChIP == "has_data", copyNumber == "sCopy") %>% 
  dplyr::arrange(SM_ID, SMTF)

productionDataList <- purrr::transpose(productionData) %>% 
  purrr::set_names(nm = purrr::map(., "SMTF"))

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% productionData$degId)


rnaseqInfoList <- purrr::transpose(rnaseqInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "comparison"))

combinedDegs <- suppressMessages(readr::read_tsv(file = file_polIIDegs))


## extract genes belonging to "GO:0022900 electron transport chain"
geneset <- AnnotationDbi::select(
  x = orgDb, keys = "GO:0022900", columns = c("GID", "GENE_NAME"), keytype = "GOALL"
) %>% 
  dplyr::select(geneId = GID, geneName = GENE_NAME)

genePos <- as.data.frame(GenomicFeatures::genes(txDb)) %>%
  dplyr::select(geneId = gene_id, chr = seqnames, start, end, strand)

geneset <- dplyr::left_join(x = geneset, y = genePos, by = "geneId") %>%
  dplyr::arrange(chr, start) %>%
  dplyr::mutate(
    geneId = forcats::as_factor(geneId)
  )

###########################################################################
## get DEG data for each SM_TF cluster from its own polII_DEG set
degData <- NULL
rowId <- 1

for (rowId in 1:nrow(productionData)) {
  
  smTf <- productionData$SMTF[rowId]
  smTfName <- productionData$SMTF_name[rowId]
  degId <- productionData$degId[rowId]
  
  degs <- dplyr::filter(combinedDegs, comparison == degId) %>% 
    dplyr::select(geneId, !!col_lfc, shrinkLog2FC, pvalue, padj) %>% 
    dplyr::mutate(
      comparison = degId,
      SMTF = smTf,
      SMTF_name = smTfName
    )
  
  degData <- dplyr::left_join(
    x = geneset, y = degs, by = "geneId"
  ) %>% 
    dplyr::bind_rows(degData)
  
}


readr::write_tsv(x = degData, file = paste(outPrefix, ".combined_DEGs.tab", sep = ""))

###########################################################################

degData <- dplyr::mutate(
  degData,
  newLfc = dplyr::if_else(
    condition = padj <= 0.05, true = !!sym(col_lfc), false = 0, missing = 0
  )
)

lfcMat <- tidyr::pivot_wider(
  data = degData,
  id_cols = c(SMTF_name),
  names_from = geneName ,
  values_from = newLfc
) %>% 
  tibble::column_to_rownames(var = "SMTF_name") %>% 
  dplyr::select(!!!geneset$geneName) %>% 
  as.matrix()


fcHeatmap <- Heatmap(
  matrix = lfcMat,
  col = colorRamp2(breaks = c(-3, 0, 3), colors = c("#313695", "gray", "#d73027"), space = "LAB"),
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  cluster_columns = FALSE,
  show_row_names = TRUE,
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



