library(tidyverse)


rm(list = ls())

###########################################################################
#######################
### DESeq2 pairwise ###
#######################

## RNAseq DESeq2 differential gene expression batches
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
diffDataPath <- here::here("analysis", "07_polII_diff")
script_deseq2 <- here::here("scripts", "13_DESeq2_diff.R")

runConfig <- suppressMessages(readr::read_tsv(file_RNAseq_info))

i <- 1

for (i in 1:nrow(runConfig)) {
  cat("Processing DESeq2:", runConfig$comparison[i], "...\n")
  
  system2(
    command = "Rscript",
    args = c(script_deseq2, "--config", file_RNAseq_info, "--deg", runConfig$comparison[i]),
    stdout = paste("logs/stdouterr.DESeq2.", runConfig$comparison[i],".log", sep = ""),
    stderr = paste("logs/stdouterr.DESeq2.", runConfig$comparison[i],".log", sep = "")
  )
  
}

###########################################################################
####################################
### DESeq2 Functional enrichment ###
####################################

file_RNAseq_info <- here::here("data", "DESeq2_DEG_info.txt")
diffDataPath <- here::here("analysis", "02_DESeq2_diff")
script_deseq2_GO <- here::here("scripts", "03_RNAseq_functional_enrichment.R")

runConfig <- suppressMessages(readr::read_tsv(file_RNAseq_info))

i <- 1

for (i in 1:nrow(runConfig)) {
  
  system2(
    command = "Rscript",
    # args = c(script_deseq2_GO, "--help"),
    args = c(script_deseq2_GO, "--config", file_RNAseq_info, "--deg", runConfig$comparison[i]),
    stdout = paste("stdouterr.functional_enrichment.", runConfig$comparison[i],".log", sep = ""),
    stderr = paste("stdouterr.functional_enrichment.", runConfig$comparison[i],".log", sep = "")
  )
  
}

###########################################################################



## copy results for sharing
file_RNAseq_info <- here::here("data", "DESeq2_DEG_info.txt")
diffDataPath <- here::here("analysis", "02_DESeq2_diff")

runConfig <- suppressMessages(readr::read_tsv(file_RNAseq_info))

tempStore <- here::here("analysis", "temp_store")
dir.create(tempStore)
i <- 1

for (i in 1:nrow(runConfig)) {
  
  newDir <- paste(tempStore, "/", runConfig$comparison[i], sep = "")
  dir.create(newDir)
  
  file_deseq2Res <- paste(diffDataPath, "/", runConfig$comparison[i], "/", runConfig$comparison[i], ".DEG_all.xlsx", sep = "")
  file_enrich <- paste(diffDataPath, "/", runConfig$comparison[i], "/", runConfig$comparison[i], ".enrichment.xlsx", sep = "")
  file_plots <- paste(diffDataPath, "/", runConfig$comparison[i], "/", runConfig$comparison[i], ".summary_plots.pdf", sep = "")

  file.copy(from = file_deseq2Res, to = newDir)
  file.copy(from = file_enrich, to = newDir)
  file.copy(from = file_plots, to = newDir)
  
}




