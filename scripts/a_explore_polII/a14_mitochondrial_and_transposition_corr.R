suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))



rm(list = ls())

source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################

analysisName <- "respiration_vs_transposition"
outDir <- here::here("analysis", "08_polII_analysis", "mitochondria")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
diffDataPath <- here::here("analysis", "06_polII_diff")
file_sampleInfo <- here::here("data", "reference_data", "polII_sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_polIIDegs <- here::here("analysis", "06_polII_diff", "polII_DEGs.fpkm_filtered.combined.tab")

file_peakDegStats <- here::here("analysis", "10_TF_polII_integration", "stats", "peak_DEG_stats.tab")

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

peakDegStats <- suppressMessages(readr::read_tsv(file_peakDegStats))

## extract genes belonging to "GO:0022900 electron transport chain"
geneset <- AnnotationDbi::select(
  x = orgDb, keys = c("GO:0022900", "GO:0015074", "GO:0006313"),
  columns = c("GID", "GENE_NAME"), keytype = "GOALL"
) %>% 
  dplyr::select(geneId = GID, geneName = GENE_NAME, GOID = GOALL)

goDesc <- AnnotationDbi::select(
  x = GO.db, keys = unique(geneset$GOID), columns = "TERM", keytype = "GOID"
)

geneset <- dplyr::left_join(x = geneset, y = goDesc, by = "GOID") %>% 
  tidyr::unite(col = "id", geneId, GOID, sep = "_", remove = FALSE)


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
  id_cols = c(id),
  names_from = SMTF_name,
  values_from = newLfc
) %>% 
  tibble::column_to_rownames(var = "id") %>% 
  as.matrix()

rowAnDf <- tibble::tibble(id = rownames(lfcMat)) %>% 
  dplyr::left_join(y = geneset, by = "id")

colAnDf <- tibble::tibble(SMTF_name = colnames(lfcMat)) %>% 
  dplyr::left_join(y = peakDegStats, by = "SMTF_name")

topAn <- ComplexHeatmap::HeatmapAnnotation(
  peaks = anno_barplot(
    x = dplyr::select(colAnDf, peaks) %>% tibble::deframe(),
    gp = gpar(fill = "#FFC20A", col = "#FFC20A"),
    height = unit(2, "cm")
  ),
  degs = anno_barplot(
    x = dplyr::select(colAnDf, up, down) %>% as.matrix(),
    gp = gpar(fill = c("#a50026", "#313695"), col = c("#a50026", "#313695")),
    height = unit(2, "cm")
  ),
  annotation_label = c("#peaks", "#DEGs"),
  annotation_name_gp = gpar(fontsize = 14, fontface = "bold")
)

fcHeatmap <- ComplexHeatmap::Heatmap(
  matrix = lfcMat,
  col = colorRamp2(breaks = c(-3, 0, 3), colors = c("#313695", "gray", "#d73027"), space = "LAB"),
  row_split = rowAnDf$TERM,
  top_annotation = topAn,
  row_title_rot = 0,
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 18, fontface = "bold"), 
  cluster_rows = TRUE, show_row_dend = FALSE,
  cluster_columns = TRUE, show_column_dend = FALSE,
  show_row_names = FALSE, show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 16),
  # width = unit(4, "cm"),
  heatmap_legend_param = list(title = "\nlog2(OE/WT)")
)


pdf(file = paste(outPrefix, ".pdf", sep = ""), width = 16, height = 10)

draw(
  object = fcHeatmap,
  column_title = "Correlation between mitochondria count and DNA transposition",
  column_title_gp = gpar(fontsize = 20, fontface = "bold")
)

dev.off()


