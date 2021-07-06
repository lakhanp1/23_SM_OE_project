suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(require(openxlsx))


## draw heatmap of summarized GO terms and their enrichment ratio in SMTFOE DEGs

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/04_GO_enrichment/s01_enrichment_functions.R")
source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################
analysisName <- "DEG_GO_summary"
outDir <- here::here("analysis", "08_polII_analysis", "04_DEG_func_enrich")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_polIIDegs <- here::here("analysis", "06_polII_diff", "polII_DEGs.fpkm_filtered.combined.tab")

file_revigoSummary <- here::here(
  "analysis", "08_polII_analysis", "04_DEG_func_enrich", "REVIGO_summarized_GO.txt"
)

TF_dataPath <- here::here("data", "TF_data")
diffDataPath <- here::here("analysis", "06_polII_diff")

orgDb <- org.Anidulans.FGSCA4.eg.db

cutoff_macs2Pval <- 20

col_lfc <- "log2FoldChange"
cutoff_fpkm <- 10
cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1
##################################################################################

productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::filter(has_polII_ChIP == "has_data", has_TF_ChIP == "has_data", copyNumber == "sCopy") %>% 
  dplyr::arrange(SM_ID, geneId)

productionData$OESMTF_name <- AnnotationDbi::mapIds(
  x = orgDb, keys = productionData$geneId, column = "GENE_NAME", keytype = "GID"
)

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% productionData$degId)


rnaseqInfoList <- purrr::transpose(rnaseqInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "comparison"))

combinedDegs <- suppressMessages(readr::read_tsv(file = file_polIIDegs))
goIds <- suppressMessages(readr::read_tsv(file = file_revigoSummary))

##################################################################################
## get topGO data for each DEG set
rowId <- 1
mergedData <- NULL

for (rowId in 1:nrow(productionData)) {
  
  tfSampleId <- productionData$tfId[rowId]
  degId <- productionData$degId[rowId]
  
  degs <- dplyr::filter(combinedDegs, comparison == degId)
  
  downDegs <- degs %>% 
    dplyr::filter(
      padj <= cutoff_fdr, !!sym(col_lfc) <= cutoff_down, fpkmFilter == "pass"
    ) %>% 
    dplyr::mutate(category = "down")
  
  upDegs <- degs %>% 
    dplyr::filter(
      padj <= cutoff_fdr & !!sym(col_lfc) >= cutoff_up, fpkmFilter == "pass"
    ) %>% 
    dplyr::mutate(category = "up")
  
  
  downGoMap <- GO_map(
    genes = downDegs$geneId, goTerms = goIds$id, orgdb = orgDb, keytype = "GID"
  ) %>% 
    dplyr::mutate(
      geneRatio = -1*enrichment,
      category = "down"
    )
  
  upGoMap <- GO_map(
    genes = upDegs$geneId, goTerms = goIds$id, orgdb = orgDb, keytype = "GID"
  ) %>% 
    dplyr::mutate(
      geneRatio = enrichment,
      category = "up"
    )
  
  degGoMap <- dplyr::bind_rows(upGoMap, downGoMap) %>% 
    dplyr::mutate(
      comparison = degId,
      OESMTF = !!productionData$geneId[rowId],
      OESMTF_name = !!productionData$OESMTF_name[rowId]
    ) %>% 
    dplyr::filter(count != 0) %>% 
    dplyr::group_by(GOID) %>% 
    dplyr::arrange(desc(count), .by_group = FALSE) %>% 
    dplyr::slice(1L) %>% 
    dplyr::ungroup()
  
  mergedData <- dplyr::bind_rows(mergedData, degGoMap)
  
}

readr::write_tsv(x = mergedData, file = paste(outPrefix, ".data2.tab", sep = ""))

wideDf <- tidyr::pivot_wider(
  data = mergedData,
  id_cols = c(GOID, TERM),
  names_from = OESMTF_name,
  values_from = geneRatio,
  values_fill = 0
)


dataMat <- dplyr::select(wideDf, -TERM) %>% 
  tibble::column_to_rownames(var = "GOID") %>% 
  as.matrix()

ht_go <- Heatmap(
  matrix = dataMat,
  col = colorRamp2(breaks = c(-0.1, 0, 0.1), colors = c("#313695", "gray", "#d73027"), space = "LAB"),
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  cluster_columns = TRUE,
  show_row_names = TRUE,
  row_labels = wideDf$TERM,
  row_names_gp = gpar(fontsize = 16),
  row_names_max_width = unit(18, "cm"),
  column_names_gp = gpar(fontsize = 16), 
  heatmap_legend_param = list(title = "\ngene ratio")
)

pdf(file = paste(outPrefix, ".heatmap.pdf", sep = ""), height = 15, width = 22)

draw(
  object = ht_go,
  column_title = "GO term summary of SMTF-OE/WT DEGs",
  column_title_gp = gpar(fontsize = 20, fontface = "bold")
)

dev.off()
