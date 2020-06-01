suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(BSgenome.Anidulans.FGSCA4.AspGD))
suppressPackageStartupMessages(library(here))


## peak annotation pie chart for AN0153 peaks

rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

##################################################################################

analysisName <- "AN0153_ChIP_summary"
outDir <- here("analysis", "04_AN0153_analysis")


outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")

TF_dataPath <- here::here("data", "TF_data")

samples_tf <- "AN0153_sCopy_OE_16h_HA_ChIPMix46_1"


if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

orgDb <- org.Anidulans.FGSCA4.eg.db
genome <- BSgenome.Anidulans.FGSCA4.AspGD
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

##################################################################################


tfInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = samples_tf,
                                   dataPath = TF_dataPath)

smGenes <- suppressMessages(
  AnnotationDbi::select(
    x = orgDb,
    keys = keys(x = orgDb, keytype = "SM_GENE"),
    columns = c("SM_GENE", "SM_CLUSTER", "SM_ID"),
    keytype = "SM_GENE")) %>% 
  dplyr::mutate(geneType = "SM")

i <- 1
peakType <- dplyr::case_when(
  tfInfo$peakType[i] == "narrow" ~ "narrowPeak",
  tfInfo$peakType[i] == "broad" ~ "broadPeak"
)

backboneGene <- tfInfo$SM_TF[i]

## few TFs are mapped to multiple SM clusters. So preparing the list of SM tf data
smTfInfo <- suppressMessages(
  AnnotationDbi::select(
    x = orgDb,
    keys = na.omit(backboneGene),
    columns = c("SM_GENE", "SM_CLUSTER", "SM_ID"),
    keytype = "GID")) %>% 
  dplyr::filter(!is.na(SM_ID))

clusterGenes <- suppressMessages(
  AnnotationDbi::select(x = orgDb,
                        keys = smTfInfo$SM_ID,
                        columns = c("GID", "SM_GENE"),
                        keytype = "SM_ID")
)

genesToMark <- list(
  SM_cluster = unique(clusterGenes$SM_GENE),
  other_SM_genes = setdiff(smGenes$SM_GENE, clusterGenes$SM_GENE)
)


chipSummary <- chip_summary(
  sampleId = tfInfo$sampleId[i],
  peakAnnotation = tfInfo$peakAnno[i],
  peakFile = tfInfo$peakFile[i],
  peakFormat = peakType,
  markTargets = genesToMark,
  pointColor = structure(c("red", "blue"), names = names(genesToMark)),
  pointAlpha = structure(c(1, 0.5), names = names(genesToMark))
)

svg(filename = paste(outPrefix, ".peak_annotation.pie.svg", sep = ""), width = 10, height = 10)
chipSummary$plots$annoPie
dev.off()









