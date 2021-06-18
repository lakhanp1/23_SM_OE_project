suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))


## post-process RNA polII-DESeq2 DEGs to include FPKM values and filter genes
## which show low FPKM values

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

##################################################################################

diffDataPath <- here::here("analysis", "06_polII_diff")
file_sampleInfo <- here::here("data", "reference_data", "polII_sample_info.txt")
file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_fpkm <- here::here("data", "polII_data", "polII_signal_matrix.tab")

orgDb <- org.Anidulans.FGSCA4.eg.db

cutoff_fpkm <- 10
cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1

##################################################################################

productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::filter(has_polII_ChIP == "has_data")

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% productionData$degId)

rnaseqInfoList <- purrr::transpose(rnaseqInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "comparison"))

fpkmMat <- suppressMessages(readr::read_tsv(file = file_fpkm))

##################################################################################

mergedData <- NULL

rowId <- 1

for (rowId in 1:nrow(rnaseqInfo)) {
  
  degInfo <- purrr::pluck(.x = rnaseqInfoList, rnaseqInfo$comparison[[rowId]])
  
  degDir <- paste(diffDataPath, "/", degInfo$comparison, sep = "")
  outPrefix <- paste(degDir, "/", degInfo$comparison, sep = "")
  sampleIds <- unlist(stringr::str_split(string = degInfo$samples, pattern = ";"))
  
  degData <- suppressMessages(
    readr::read_tsv(
      file = paste(degDir, "/", degInfo$comparison, ".DEG_all.txt", sep = "")
    )
  )
  
  fpkmData <- dplyr::select(fpkmMat, geneId, !!!sampleIds) %>% 
    dplyr::mutate(
      maxFpkm = pmax(!!!syms(sampleIds)),
      minFpkm = pmin(!!!syms(sampleIds)),
      sampleIds = paste(!!!sampleIds, sep = ";")
    ) %>% 
    tidyr::unite(col = fpkms, !!!syms(sampleIds), sep = ";")
  
  degData <- dplyr::left_join(
    x = degData, y = fpkmData, by = "geneId"
  ) %>% 
    dplyr::mutate(
      fpkmFilter = dplyr::if_else(
        condition = maxFpkm >= cutoff_fpkm, true = "pass", false = "fail", missing = "fail"
      )
    )
  
  mergedData <- dplyr::bind_rows(mergedData, degData)
  ########################
  ## store data
  readr::write_tsv(
    x = degData, file = paste(outPrefix, ".DEG_FPKM_filtered.txt", sep = "")
  )
  
  ## write data to excel file
  wb <- openxlsx::createWorkbook(creator = "Lakhansing Pardehi")
  openxlsx::addWorksheet(wb = wb, sheetName = "DEG_FPKM_filtered")
  openxlsx::writeData(
    wb = wb, sheet = 1, startCol = 2, startRow = 1,
    x = paste("## DESeq2 results with FPKM filtering information:", degInfo$comparison)
  )
  openxlsx::writeData(
    wb = wb, sheet = 1, x = degData,
    startCol = 1, startRow = 2, withFilter = TRUE,
    keepNA = TRUE, na.string = "NA"
  )
  headerStyle <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#e6e6e6")
  openxlsx::addStyle(wb = wb, sheet = 1, style = headerStyle, rows = 2, cols = 1:ncol(degData))
  openxlsx::setColWidths(wb = wb, sheet = 1, cols = 1, widths = "auto")
  openxlsx::freezePane(wb = wb, sheet = 1, firstActiveRow = 3, firstActiveCol = 2)
  
  # openxlsx::openXL(wb)
  openxlsx::saveWorkbook(
    wb = wb,
    file = paste(outPrefix, ".DEG_FPKM_filtered.xlsx", sep = ""),
    overwrite = TRUE
  )
  
}


readr::write_tsv(
  x = mergedData,
  file = paste(diffDataPath, "/polII_DEGs.fpkm_filtered.combined.tab", sep = "")
)


