library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(TxDb.Anidulans.FGSCA4.AspGD.GFF)
library(here)
library(ggbeeswarm)
library(ggpubr)
library(ggrepel)


rm(list = ls())


##################################################################################
analysisName <- "PolII_replicate"
outDir <- here::here("analysis", "02_QC_polII")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_replicates <- here::here("analysis", "02_QC_polII", "polII_replicates.txt")

file_exptInfo <- here::here("data", "referenceData/sample_info.txt")

file_genes <- here::here("data", "referenceData/AN_genes_for_polII.bed")
orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

polII_dataPath <- here::here("data", "polII_data")

file_polII <- paste(polII_dataPath, "/", "sample_polII.list", sep = "")

##################################################################################

geneSet <- suppressMessages(readr::read_tsv(
  file = file_genes,
  col_names = c("chr", "start", "end", "geneId", "score", "strand"))) %>%
  dplyr::select(geneId)

polIIsamples <- readr::read_tsv(file = file_polII, col_names = c("sampleId"),  comment = "#")

polIIInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = polIIsamples$sampleId,
  dataPath = polII_dataPath)


replicatePairs <- suppressMessages(readr::read_tsv(file = file_replicates))

plotListAll <- list()
plotList_pval_distibution <- list()

i <- 1

pdf(file = paste(outPrefix, ".correlation.pdf", sep = ""), width = 15, height = 10,
    onefile = TRUE, pointsize = 8)


for (i in 1:nrow(replicatePairs)) {
  cat(i, ":", replicatePairs$rep1[i], "\n")
  
  repInfo <- dplyr::filter(polIIInfo, sampleId %in% c(replicatePairs$rep1[i], replicatePairs$rep2[i]))
  
  if(nrow(repInfo) == 0){
    next
  }
  
  if(nrow(repInfo) != 2){
    stop("sampleInfo should have 2 rows for 2 replicates to be compared")
  }
  
  rep1Col <- repInfo$sampleId[1]
  rep2Col <- repInfo$sampleId[2]
  exprDf <- get_polII_expressions(genesDf = geneSet, exptInfo = repInfo) 
  
  pairCor <- compare_replicates(data = exprDf, rep1Col = rep1Col, rep2Col = rep2Col,
                                trans = "log2")
  
  plot(pairCor$figure)
  
}

dev.off()

##################################################################################









