suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(openxlsx))



## process RankComp2 output and prepare result table

rm(list = ls())

##################################################################################

rankcompPath <- here::here("analysis", "07_polII_rank_diff")

file_sampleInfo <- here::here("data", "reference_data", "polII_sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_fpkm <- here::here("data", "polII_data", "polII_signal_matrix.tab")
file_geneIndex <- paste(rankcompPath, "/gene_index.tab", sep = "")

cutoff_rank <- 5
cutoff_rankUp <- cutoff_rank
cutoff_rankDown <- cutoff_rank * -1

##################################################################################

diffInfo <- suppressMessages(readr::read_tsv(file = file_RNAseq_info))

sampleInfo <- read.table(file = file_sampleInfo, header = T, sep = "\t", stringsAsFactors = F)
fpkmMat <- suppressMessages(readr::read_tsv(file = file_fpkm))

geneIndex <- suppressMessages(readr::read_tsv(file = file_geneIndex, col_names = "geneId"))

##################################################################################
nGenes <- nrow(fpkmMat)

i <- 1

for (i in 1:nrow(diffInfo)) {
  
  contrastId <- diffInfo$comparison[i]
  
  ########################
  ## process RankComp data
  degDir <- paste(rankcompPath, "/", contrastId, sep = "")
  outPrefix <- paste(degDir, "/", contrastId, sep = "")
  
  upDeg <- suppressMessages(
    readr::read_tsv(file = paste(degDir, "/up_regulated_1_0.dat", sep = ""),
                    col_names = c("index", "fdr"))
  ) %>% 
    dplyr::mutate(diff = "up")
  
  downDeg <- suppressMessages(
    readr::read_tsv(file = paste(degDir, "/down_regulated_1_0.dat", sep = ""),
                    col_names = c("index", "fdr"))
  ) %>% 
    dplyr::mutate(diff = "down")
  
  degs <- dplyr::bind_rows(upDeg, downDeg) %>% 
    dplyr::mutate(index = index + 1)
  
  
  diffData <- geneIndex %>% 
    dplyr::mutate(contrast = contrastId, fdr = NA, diff = "noDEG")
  
  diffData$fdr[degs$index] <- degs$fdr
  diffData$diff[degs$index] <- degs$diff
  
  ## generate rank data for the samples being compared
  samples <- unlist(stringr::str_split(string = diffInfo$samples[i], pattern = ";"))
  
  mtSampleInfo <- dplyr::filter(sampleInfo, !!sym(diffInfo$design[i]) %in% diffInfo$group1[i]) %>% 
    dplyr::filter(sampleId %in% samples) %>% 
    dplyr::mutate(index = paste("MT", row_number(), sep = ""))
  
  wtSampleInfo <- dplyr::filter(sampleInfo, !!sym(diffInfo$design[i]) %in% diffInfo$group2[i]) %>% 
    dplyr::filter(sampleId %in% samples) %>% 
    dplyr::mutate(index = paste("WT", row_number(), sep = ""))
  
  degsetList <- dplyr::bind_rows(mtSampleInfo, wtSampleInfo) %>% 
    dplyr::mutate(rankCol = paste("rank.", sampleId, sep = "")) %>% 
    purrr::transpose() %>% 
    purrr::set_names(nm = purrr::map(., "index"))
  
  ########################
  ## build the FPKM score rank dataframe
  ranksDf <- purrr::map_dfc(
    .x = degsetList,
    .f = function(x){
      list(rank(fpkmMat[[x$sampleId]])) %>%
        purrr::set_names(nm = x$rankCol)
    }
  ) %>% 
    dplyr::mutate(geneId = fpkmMat$geneId) %>% 
    dplyr::select(geneId, everything())
  
  if(!all(ranksDf[[degsetList$WT1$rankCol]] == rank(fpkmMat[[degsetList$WT1$sampleId]]))){
    stop("ranks in does not match...")
  }
  
  
  diffPairs <- tidyr::crossing(wt = wtSampleInfo$index, mt = mtSampleInfo$index) %>% 
    tidyr::unite(col = "cmp", mt, wt, sep = "_vs_", remove = FALSE)
  
  ## compare ranks (as % difference) for different combinations of MT and WT samples
  rankCompDf <- purrr::transpose(diffPairs) %>% 
    purrr::map_dfc(
      .f = function(x){
        mtRanks <- dplyr::pull(.data = ranksDf, var = degsetList[[x$mt]]$rankCol)
        wtRanks <- dplyr::pull(.data = ranksDf, var = degsetList[[x$wt]]$rankCol)
        
        list(round((mtRanks - wtRanks)*100/nGenes, digits = 3)) %>% 
          purrr::set_names(nm = x$cmp)
      }
    ) %>% 
    dplyr::mutate(
      geneId = fpkmMat$geneId,
      minRankDiff = pmin(!!!syms(diffPairs$cmp)),
      maxRankDiff = pmax(!!!syms(diffPairs$cmp))
    ) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(
      avgRankDiff = round(mean(c_across(contains("_vs_"))), digits = 3),
      sdRankDiff = round(sd(c_across(contains("_vs_"))), digits = 3),
      consensus = dplyr::case_when(
        all(sign(c_across(contains("_vs_"))) == 1) ~ "+",
        all(sign(c_across(contains("_vs_"))) == -1) ~ "-",
        TRUE ~ "+-"
      ),
      diffStrong = dplyr::case_when(
        consensus == "-" & maxRankDiff < cutoff_rankDown ~ "down",
        consensus == "-" & minRankDiff > cutoff_rankUp ~ "up",
        
      )
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(geneId, everything())
  
  
  rankCompData <- dplyr::left_join(x = diffData, y = ranksDf, by = "geneId") %>% 
    dplyr::left_join(y = rankCompDf, by = "geneId")
  
  ########################
  ## store data
  write_tsv(x = rankCompData, path = paste(outPrefix, ".RankComp2.tab", sep = ""))
  
  ## write data to excel file
  wb <- openxlsx::createWorkbook(creator = "Lakhansing Pardehi Genomics Core")
  openxlsx::addWorksheet(wb = wb, sheetName = "RankComp2")
  openxlsx::writeData(
    wb = wb, sheet = 1, startCol = 2, startRow = 1,
    x = paste("## Differential gene expression rank analysis by RankComp2 for",
              contrastId,":", diffInfo$group1[i], "/", diffInfo$group2[i])
  )
  openxlsx::writeData(
    wb = wb, sheet = 1, x = rankCompData,
    startCol = 1, startRow = 2, withFilter = TRUE,
    keepNA = TRUE, na.string = "NA"
  )
  headerStyle <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#e6e6e6")
  openxlsx::addStyle(wb = wb, sheet = 1, style = headerStyle, rows = 2, cols = 1:ncol(rankCompData))
  openxlsx::setColWidths(wb = wb, sheet = 1, cols = 1, widths = "auto")
  openxlsx::setColWidths(wb = wb, sheet = 1, cols = 2, widths = str_length(contrastId))
  openxlsx::freezePane(wb = wb, sheet = 1, firstActiveRow = 3, firstActiveCol = 2)
  
  # openxlsx::openXL(wb)
  openxlsx::saveWorkbook(wb = wb, file = paste(outPrefix, ".RankComp2.xlsx", sep = ""), overwrite = TRUE)
  
}












