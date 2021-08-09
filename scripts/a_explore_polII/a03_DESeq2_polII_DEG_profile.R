suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/04_GO_enrichment/s01_enrichment_functions.R")
source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

##################################################################################
analysisName <- "polII_DEG_profile_heatmaps"
outDir <- here::here("analysis", "08_polII_analysis", "05_DEG_profile_heatmap")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_polIIDegs <- here::here("analysis", "06_polII_diff", "polII_DEGs.fpkm_filtered.combined.tab")

polII_dataPath <- here::here("data", "polII_data")
diffDataPath <- here::here("analysis", "06_polII_diff")

cutoff_fpkm <- 10
cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- -1 * cutoff_lfc

col_lfc <- "log2FoldChange"

orgDb <- org.Anidulans.FGSCA4.eg.db

matrixType <- "normalizedmatrix"
matrixDim = c(200, 200, 100, 10)

###########################################################################

productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::filter(has_polII_ChIP == "has_data", has_TF_ChIP == "has_data", copyNumber == "sCopy") %>% 
  dplyr::filter(!is.na(SM_ID)) %>% 
  dplyr::arrange(SM_ID, SMTF)

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% productionData$degId)

rnaseqInfoList <- purrr::transpose(rnaseqInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "comparison"))

combinedDegs <- suppressMessages(readr::read_tsv(file = file_polIIDegs))

###########################################################################

rowId <- 3

pdf(file = paste(outPrefix, ".pdf", sep = ""), width = 10, height = 8, onefile = TRUE)

for (rowId in 1:nrow(productionData)) {
  
  degSubInfo <- purrr::pluck(rnaseqInfoList, productionData$degId[rowId])
  
  degs <- dplyr::filter(combinedDegs, comparison == degSubInfo$comparison) %>% 
    dplyr::mutate(
      category = dplyr::case_when(
        padj <= cutoff_fdr & !!sym(col_lfc) <= cutoff_down ~ "down",
        padj <= cutoff_fdr & !!sym(col_lfc) >= cutoff_up ~ "up",
        TRUE ~ "noDEG"
      ),
      rankMetric = (-log10(pvalue) * sign(shrinkLog2FC))
    ) %>% 
    dplyr::filter(
      padj <= cutoff_fdr, abs(!!sym(col_lfc)) > cutoff_lfc, fpkmFilter == "pass"
    ) %>% 
    dplyr::mutate(
      category = forcats::fct_relevel(category, "up", "down")
    )
  
  sampleIds <- unlist(stringr::str_split(string = degSubInfo$samples, pattern = ";"))
  
  polII_info <- get_sample_information(
    exptInfoFile = file_exptInfo,
    samples = sampleIds,
    dataPath = polII_dataPath) %>% 
    dplyr::mutate(
      sampleName = paste(gene, "_polII_", rep, sep = "")
    )
  
  
  matList <- chipmine::import_profiles(
    exptInfo = polII_info,
    source = matrixType,
    up = matrixDim[1], target = matrixDim[2], down = matrixDim[3]
  )
  
  
  ## polII colors
  polIIMeanProfile <- NULL
  polIIColorList <- NULL
  if(nrow(polII_info) == 1){
    polIIMeanProfile <- matList[[polII_info$sampleId]]
  } else{
    polIIMeanProfile <- getSignalsFromList(lt = matList[polII_info$sampleId])
  }
  quantile(polIIMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
  polIIMeanColor <- colorRamp2(
    breaks = quantile(polIIMeanProfile, c(0.01, 0.5, 0.995), na.rm = T),
    colors = c("blue", "white", "red")
  )
  polIIColorList <- sapply(X = polII_info$sampleId, FUN = function(x){return(polIIMeanColor)})
  
  colorList <- unlist(list(polIIColorList))
  
  clusterColors = c("down" = "#313695", "noDEG" = "#F7F7F7", "up" = "#d73027")
  
  # ylimList <- list()
  # ylimList <- append(
  #   x = ylimList,
  #   values = sapply(c(tfIds, inputIds), function(x){return(c(0, 25))}, simplify = FALSE)
  # )
  
  multiProfiles <- chipmine::multi_profile_plots(
    exptInfo = polII_info,
    genesToPlot = degs$geneId,
    matSource = matrixType,
    matBins = matrixDim,
    # cluster_rows = TRUE, clustering_distance_rows = dist_by_closeness,
    clusters = dplyr::select(degs, geneId, cluster = category),
    profileColors = colorList, clusterColor = clusterColors,
    column_title_gp = gpar(fontsize = 12)
  )
  
  deg_htlist <- multiProfiles$heatmapList
  
  # rowOrd <- order(degs$log2FoldChange, decreasing = TRUE)
  rowOrd <- order(degs$rankMetric, decreasing = TRUE)
  
  ptTitle <- paste(productionData$SMTF_name[rowId], "-OE/WT: polII-ChIPseq DEG profile heatmap", sep = "")
  
  draw(deg_htlist,
       main_heatmap = polII_info$profileName[1],
       # annotation_legend_list = list(profile1$legend),
       column_title = ptTitle,
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       row_sub_title_side = "left",
       heatmap_legend_side = "right",
       gap = unit(7, "mm"),
       row_order = rowOrd,
       padding = unit(rep(0.5, times = 4), "cm")
  )
  
}

dev.off()


###########################################################################






