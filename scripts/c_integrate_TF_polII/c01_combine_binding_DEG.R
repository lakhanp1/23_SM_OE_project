suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))


rm(list = ls())

source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

##################################################################################

analysisName <- "binding_DEG_data"
outDir <- here::here("analysis", "10_TF_polII_integration")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_polIIDegs <- here::here("analysis", "06_polII_diff", "polII_DEGs.fpkm_filtered.combined.tab")

TF_dataPath <- here::here("data", "TF_data")
diffDataPath <- here::here("analysis", "06_polII_diff")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

cutoff_macs2Pval <- 20

cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1

col_lfc <- "log2FoldChange"
col_pval <- "pvalue"

##################################################################################

productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::filter(has_polII_ChIP == "has_data", has_TF_ChIP == "has_data", copyNumber == "sCopy") %>% 
  dplyr::arrange(SM_ID, SMTF)

productionDataList <- purrr::transpose(productionData) %>% 
  purrr::set_names(nm = purrr::map(., "SMTF"))

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = productionData$tfId,
  dataPath = TF_dataPath)

tfInfoList <- purrr::transpose(tfInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% productionData$degId)


rnaseqInfoList <- purrr::transpose(rnaseqInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "comparison"))

combinedDegs <- suppressMessages(readr::read_tsv(file = file_polIIDegs))


##################################################################################

mergedData <- NULL

rowId <- 1

for (rowId in 1:nrow(productionData)) {
  
  smTf <- productionData$SMTF[rowId]
  smTfName <- productionData$SMTF_name[rowId]
  tfId <- productionData$tfId[rowId]
  degId <- productionData$degId[rowId]
  
  cat(rowId, ": ", smTf, "\n", sep = "")
  
  degs <- dplyr::filter(combinedDegs, comparison == degId)
  
  # ## extract peak annotation
  # peakAn <- suppressMessages(readr::read_tsv(tfInfoList[[tfId]]$peakAnno)) %>% 
  #   dplyr::filter(peakPval >= cutoff_macs2Pval) %>% 
  #   dplyr::mutate(hasPeak = TRUE) %>% 
  #   dplyr::group_by(geneId) %>% 
  #   dplyr::arrange(desc(peakPval), .by_group = TRUE) %>% 
  #   dplyr::slice(1L) %>% 
  #   dplyr::ungroup() %>% 
  #   dplyr::select(geneId, peakPval, peakAnnotation, peakCategory, peakPosition, hasPeak)
  
  ## prepare TF target gene list
  peakAn <- suppressMessages(readr::read_tsv(file = tfInfoList[[tfId]]$peakAnno)) %>% 
    dplyr::filter(peakPval >= cutoff_macs2Pval) %>%
    dplyr::filter(!peakCategory %in% c("intergenic", "blacklist")) %>% 
    dplyr::filter(peakChr != "mito_A_nidulans_FGSC_A4") %>% 
    dplyr::select(
      geneId, peakId, peakRegion, peakPval, relativeSummitPos, peakAnnotation,
      peakCategory, peakPosition, relativePeakPos, summitDist, peakDist,
      bidirectional, targetOverlap, peakOverlap
    )
  
  ## combine rows when a gene has multiple peaks' annotation
  if(nrow(peakAn) > 0){
    peakAn <- dplyr::mutate(
      peakAn,
      peakCategory = forcats::fct_relevel(
        .f = peakCategory, "nearStart", "peakInFeature", "featureInPeak",
        "nearEnd", "upstreamTss"
      )
    ) %>% 
      dplyr::group_by(geneId) %>% 
      dplyr::arrange(peakCategory, .by_group = TRUE) %>% 
      dplyr::mutate(annotationRank = 1:n()) %>% 
      dplyr::summarise_at(.vars = vars(-group_cols()), .funs = ~ paste(., collapse = ";"))
    
  } else{
    ## use only geneId column to avoid left_join() error
    peakAn <- dplyr::select(peakAn, geneId)
  }
  
  
  
  bindingDegData <- dplyr::left_join(
    x = degs, y = peakAn, by = "geneId"
  ) %>% 
    dplyr::mutate(
      OESMTF = smTf,
      OESMTF_name = smTfName
    ) %>% 
    dplyr::select(
      OESMTF, OESMTF_name, geneId, geneName = GENE_NAME, everything()
    )
  
  mergedData <- dplyr::bind_rows(mergedData, bindingDegData)
  
  
}

## add filters and additional annotations  for the genes
mergedData2 <- dplyr::mutate(
  mergedData,
  significance = dplyr::case_when(
    !!sym(col_pval) <= cutoff_fdr & abs(!!sym(col_lfc)) >= cutoff_lfc & fpkmFilter == "pass" ~ "significant",
    TRUE ~ "non-significant"
  ),
  binding = dplyr::if_else(
    condition = is.na(peakId), true = "not-bound", false = "bound", missing = "not-bound"
  )
) %>% 
  dplyr::mutate(
    regulationType = dplyr::case_when(
      binding == "bound" & log2FoldChange < 0 & significance == "significant" ~ "bound+down",
      binding == "bound" & log2FoldChange > 0 & significance == "significant" ~ "bound+up",
      binding == "bound" ~ "bound",
      log2FoldChange < 0 & significance == "significant" ~ "down",
      log2FoldChange > 0 & significance == "significant" ~ "up",
      TRUE ~ "no-regulation"
    )
  )


readr::write_tsv(x = mergedData2, file = paste(outPrefix, ".merged.tab", sep = ""))

# inDf <- suppressMessages(readr::read_tsv(file = paste(outPrefix, ".merged.tab", sep = "")))




