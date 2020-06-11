suppressPackageStartupMessages(library(tidyverse))



## prepare polII FPKM data files for input to RankComp tool

rm(list = ls())

##################################################################################

outDir <- here::here("analysis", "07_polII_rank_diff")

file_sampleInfo <- here::here("data", "reference_data", "polII_sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_fpkm <- here::here("data", "polII_data", "polII_signal_matrix.tab")

##################################################################################


diffInfo <- suppressMessages(readr::read_tsv(file = file_RNAseq_info))

sampleInfo <- read.table(file = file_sampleInfo, header = T, sep = "\t", stringsAsFactors = F)
fpkmMat <- suppressMessages(readr::read_tsv(file = file_fpkm))

projDir <- getwd()

readr::write_lines(x = fpkmMat$geneId, path = paste(outDir, "/", "gene_index.tab", sep = ""))

i <- 1

## write the FPKM values in tabular format for WT and OE strain samples
for (i in 1:nrow(diffInfo)) {
  
  samples <- unlist(stringr::str_split(string = diffInfo$samples[i], pattern = ";"))
  
  oeSampleInfo <- dplyr::filter(sampleInfo, !!sym(diffInfo$design[i]) %in% diffInfo$group1[i]) %>% 
    dplyr::filter(sampleId %in% samples)
  
  wtSampleInfo <- dplyr::filter(sampleInfo, !!sym(diffInfo$design[i]) %in% diffInfo$group2[i]) %>% 
    dplyr::filter(sampleId %in% samples)
  
  diffDir <- paste(outDir, "/", diffInfo$comparison[i], sep = "")
  if(!dir.exists(diffDir)){
    dir.create(path = diffDir)
    
  }
  
  file_case <- paste(diffDir, "/", diffInfo$group1[i], ".fpkm.tab", sep = "")
  readr::write_tsv(
    x = fpkmMat[, oeSampleInfo$sampleId], col_names = FALSE,
    path = file_case
  )
  
  file_control <- paste(diffDir, "/", diffInfo$group2[i], ".fpkm.tab", sep = "")
  readr::write_tsv(
    x = fpkmMat[, wtSampleInfo$sampleId], col_names = FALSE,
    path = file_control
  )
  
  command_reoa <- paste(
    "cellcomp -v -f 0.05",
    paste(diffInfo$group2[i], ".fpkm.tab", sep = ""),
    paste(diffInfo$group1[i], ".fpkm.tab", sep = "")
  )
  
  readr::write_lines(x = command_reoa, path = paste(diffDir, "/command_reoa.sh", sep = ""))
  
}













