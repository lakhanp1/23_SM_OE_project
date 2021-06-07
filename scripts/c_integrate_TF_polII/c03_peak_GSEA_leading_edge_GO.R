suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(fgsea))


## SMTF OE peak target gene's GSEA on the SMTF OE polII DEG gene list

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")
source("D:/work_lakhan/github/omics_utils/04_GO_enrichment/s01_enrichment_functions.R")

##################################################################################

analysisName <- "DEG_peak_GSEA.leading_edge.GO"
outDir <- here::here("analysis", "10_TF_polII_integration", "peakset_enrichmet_in_DEG")
outPrefix <- paste(outDir, "/", analysisName, sep = "")
gseaPtOutDir <- paste(outDir, "/GSEA_plots", sep = "")
file_topGO <- "E:/Chris_UM/Database/A_Nidulans/annotation_resources/geneid2go.ANidulans.topGO.map"

file_gseaRes <- paste(outDir, "/peakset_enrichmet_in_DEG.fgsea.tab", sep = "")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

keggOrg <- 'ani'									                  ## KEGG organism code
col_degIdKeytype <- "GID"               ## org.db keytype of DEG file geneId
col_kegg <- "NCBI_ID"																## org.db column for NCBI ID
col_gsea <- "GID"																## org.db column to use for gsea analysis
col_topGO <- "GID"															## org.db keytype of topGO map file geneId
col_geneName <- "GENE_NAME"													## org.db column for Gene name

##################################################################################

smTf <- "AN7896"

gseaRes <- suppressMessages(readr::read_tsv(file = file_gseaRes)) %>% 
  dplyr::mutate(
    leadingEdge = stringr::str_split(string = leadingEdge, pattern = ";")
  ) %>% 
  dplyr::filter(pathway == smTf)




## topGO GO enrichment
leadingEdgeGO <- topGO_enrichment(
  goMapFile = file_topGO,
  genes = as.character(unlist(gseaRes$leadingEdge)),
  type = "BP", goNodeSize = 10,
  orgdb = orgDb, keytype = col_degIdKeytype,
  topgoColumn = col_topGO, geneNameColumn = col_geneName
)

readr::write_tsv(x = leadingEdgeGO, file = paste(outPrefix, ".", smTf, ".tab", sep = ""))

plotTitle <- paste(
  "Leading edge genes for SMTF ", smTf, " in OE/WT DEGs ", "(n = ",
  leadingEdgeGO$inputSize[1], ")", sep = ""
)

pt_leadingEdgeGO <- enrichment_scatter(df = leadingEdgeGO[1:30,], title = plotTitle)

ggsave(filename = paste(outPrefix, ".", smTf, ".png", sep = ""), plot = pt_leadingEdgeGO, width = 10, height = 10)

