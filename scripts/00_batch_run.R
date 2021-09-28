library(tidyverse)

rm(list = ls())


###########################################################################
#######################
### DESeq2 pairwise ###
#######################

# ## RNAseq DESeq2 differential gene expression batches
# file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
# diffDataPath <- here::here("analysis", "06_polII_diff")
# script_deseq2 <- here::here("scripts", "13_DESeq2_diff.R")
# 
# runConfig <- suppressMessages(readr::read_tsv(file_RNAseq_info))
# 
# rowId <- 45
# 
# for (rowId in 45:nrow(runConfig)) {
#   cat("Processing DESeq2:", runConfig$comparison[rowId], "...\n")
# 
#   system2(
#     command = "Rscript",
#     args = c(script_deseq2, "--config", file_RNAseq_info, "--deg", runConfig$comparison[rowId]),
#     stdout = paste("logs/stdouterr.DESeq2.", runConfig$comparison[rowId],".log", sep = ""),
#     stderr = paste("logs/stdouterr.DESeq2.", runConfig$comparison[rowId],".log", sep = "")
#   )
# 
# }

###########################################################################
####################################
### DESeq2 Functional enrichment ###
####################################

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
diffDataPath <- here::here("analysis", "06_polII_diff")
script_deseq2_GO <- here::here("scripts", "a_explore_polII", "a09_polII_DEGs_functional_enrichment.R")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")


productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::filter(has_polII_ChIP == "has_data")


rowId <- 1

for (rowId in 1:nrow(productionData)) {
  
  system2(
    command = "Rscript",
    # args = c(script_deseq2_GO, "--help"),
    args = c(script_deseq2_GO, "--config", file_RNAseq_info, "--deg", productionData$degId[rowId]),
    stdout = paste("logs/stdouterr.functional_enrichment.", productionData$degId[rowId],".log", sep = ""),
    stderr = paste("logs/stdouterr.functional_enrichment.", productionData$degId[rowId],".log", sep = "")
  )
  
  print(productionData$degId[rowId])
}


###########################################################################

# ## run DESeq2 and RankComp2 comparison script for all OE/WT polII comparisons
# file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
# runConfig <- suppressMessages(readr::read_tsv(file_RNAseq_info))
# script <- here::here("scripts", "21_DESeq2_RankComp.comparison.R")
# 
# rowId <- 47
# 
# for (rowId in 1:nrow(runConfig)) {
#   
#   cat("Running", basename(script), "for", runConfig$comparison[rowId], "\n")
#   
#   system2(
#     command = "Rscript",
#     # args = c(script, "--help"),
#     args = c(script, "--deg", runConfig$comparison[rowId]),
#     stdout = paste("logs/stdouterr.DESeq2_RankComp2_comp.", runConfig$comparison[rowId],".log", sep = ""),
#     stderr = paste("logs/stdouterr.DESeq2_RankComp2_comp.", runConfig$comparison[rowId],".log", sep = "")
#   )
# 
# }

###########################################################################



