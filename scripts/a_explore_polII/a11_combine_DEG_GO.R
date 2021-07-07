suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(require(openxlsx))


## combine GO enrichment results 

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/04_GO_enrichment/s01_enrichment_functions.R")
source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################
analysisName <- "DEG_GO_combined"
outDir <- here::here("analysis", "08_polII_analysis", "04_DEG_func_enrich")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")

TF_dataPath <- here::here("data", "TF_data")
diffDataPath <- here::here("analysis", "06_polII_diff")

orgDb <- org.Anidulans.FGSCA4.eg.db

##################################################################################
if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::filter(has_polII_ChIP == "has_data", has_TF_ChIP == "has_data", copyNumber == "sCopy") %>% 
  dplyr::arrange(SM_ID, SMTF)

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% productionData$degId)


rnaseqInfoList <- purrr::transpose(rnaseqInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "comparison"))

## get topGO data for each DEG set
rowId <- 1
mergedData <- NULL


for (rowId in 1:nrow(productionData)) {
  
  tfSampleId <- productionData$tfId[rowId]
  degId <- productionData$degId[rowId]
  
  goData <- suppressMessages(readr::read_tsv(file = rnaseqInfoList[[degId]]$topGO)) %>% 
    dplyr::mutate(
      OESMTF = !!productionData$SMTF[rowId],
      OESMTF_name = !!productionData$SMTF_name[rowId]
    ) %>% 
    dplyr::select(
      OESMTF, OESMTF_name, contrast, everything()
    )
  
  mergedData <- dplyr::bind_rows(mergedData, goData)
}


readr::write_tsv(
  file = paste(outPrefix, ".data.tab", sep = ""), x = mergedData
)





