suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(gginnards))

## script to test the circular visualization of SM cluster genes

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

##################################################################################

analysisName <- "SM_clusters_binding_DEG.circos"
outDir <- here::here("analysis", "10_TF_polII_integration", "03_SM_cluster_level_analysis")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_polIIDegs <- here::here("analysis", "06_polII_diff", "polII_DEGs.fpkm_filtered.combined.tab")
file_smRegions <- here::here("data", "reference_data", "SM_genes_regions.tab")

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

## extract all SM cluster genes
smGeneInfo <- AnnotationDbi::select(
  x = orgDb, keys = keys(x = orgDb, keytype = "SM_ID"),
  columns = c("GID", "GENE_NAME", "SM_CLUSTER", "TF_GENE"), keytype = "SM_ID"
) %>% 
  dplyr::rename(geneId = GID) %>% 
  dplyr::mutate(
    geneLabel = paste(geneId, "(", GENE_NAME, ")", sep = ""),
    geneLabel = if_else(condition = geneId == GENE_NAME, true = geneId, false = geneLabel)
  )

genePos <- as.data.frame(
  GenomicFeatures::genes(x = txDb, filter = list(gene_id = unique(smGeneInfo$geneId)))
) %>% 
  dplyr::select(geneId = gene_id, chr = seqnames, start, end, strand)

smGeneInfo <- dplyr::left_join(x = smGeneInfo, y = genePos, by = "geneId") %>% 
  dplyr::arrange(chr, start, SM_ID) %>% 
  dplyr::mutate(
    SM_ID = forcats::as_factor(SM_ID)
  )

smGeneRegions <- suppressMessages(readr::read_tsv(file = file_smRegions)) %>% 
  dplyr::left_join(y = genePos, by = "geneId") %>%
  dplyr::arrange(chr, start, SM_IDs)

# ##################################################################################
# ## get DEG+binding data for each SM cluster genes from all SMTFOE experiments
# rowId <- 1
# mergedData <- NULL
# 
# 
# for (rowId in 1:nrow(productionData)) {
# 
#   tfSampleId <- productionData$tfId[rowId]
#   degId <- productionData$degId[rowId]
# 
#   ## extract peak annotation
#   peakAn <- suppressMessages(readr::read_tsv(tfInfoList[[tfSampleId]]$peakAnno)) %>%
#     dplyr::filter(peakPval >= cutoff_macs2Pval) %>%
#     dplyr::mutate(hasPeak = TRUE) %>%
#     dplyr::group_by(geneId) %>%
#     dplyr::arrange(desc(peakPval), .by_group = TRUE) %>%
#     dplyr::slice(1L) %>%
#     dplyr::ungroup() %>%
#     dplyr::select(geneId, peakPval, peakAnnotation, peakCategory, peakPosition, hasPeak)
# 
#   ## extract DEG data
#   diffData <- dplyr::filter(combinedDegs, comparison == degId) %>%
#     dplyr::select(geneId, comparison, !!col_lfc, shrinkLog2FC, pvalue, padj, maxFpkm, fpkmFilter)
# 
#   bindingDegData <- dplyr::left_join(
#     x = genePos, y = diffData, by = "geneId"
#   ) %>%
#     dplyr::left_join(y = peakAn, by = "geneId") %>%
#     dplyr::mutate(
#       OESMTF = !!productionData$SMTF[rowId],
#       OESMTF_name = !!productionData$SMTF_name[rowId]
#     ) %>%
#     dplyr::select(
#       OESMTF,OESMTF_name, geneId, everything()
#     )
# 
#   mergedData <- dplyr::bind_rows(mergedData, bindingDegData)
# 
# }
# 
# 
# mergedData2 <- dplyr::mutate(
#   mergedData,
#   significance = dplyr::case_when(
#     !!sym(col_pval) <= cutoff_fdr & fpkmFilter == "pass" ~ "significant",
#     TRUE ~ "non-significant"
#   ),
#   binding = dplyr::if_else(
#     condition = hasPeak, true = "bound", false = "not-bound", missing = "not-bound"
#   ),
#   tfGene = if_else(
#     condition = !is.na(TF_GENE), true = "TF", false = "non-TF"
#   )
# )
# 
# ###########################################################################

## 310 x 220
baseData <- dplyr::arrange(.data = smGeneRegions, chr, start) %>%
  dplyr::mutate_at(
    .vars = c("chr", "geneId", "regionId"), .funs = ~forcats::as_factor(.)
  ) %>%
  dplyr::mutate(
    chr = forcats::fct_drop(chr)
  ) %>%
  dplyr::group_by(regionId) %>%
  dplyr::arrange(start, .by_group = TRUE) %>%
  dplyr::mutate(
    genePos = 1:n()
  ) %>%
  dplyr::ungroup()

pdf(file = paste(outPrefix, ".temp.pdf", sep = ""), width = 20, 20)

circos.par(gap.degree = 0, cell.padding = c(0,0,0,0))

circos.initialize(
  sectors = baseData$regionId,
  # x = baseData$genePos
  xlim = dplyr::group_by(baseData, regionId) %>% 
    dplyr::summarise(
      start = 0,
      end = n()
    ) %>% 
    tibble::column_to_rownames(var = "regionId") %>% 
    as.matrix()
)

df <- dplyr::filter(baseData, genePos %in% c(20:40)) %>%
  dplyr::mutate(y = 3)


circos.track(
  sectors = baseData$regionId,
  track.index = 1,
  x = baseData$genePos-0.5, y = rep(0, nrow(baseData)),
  ylim = c(0, 1),
  panel.fun = function(x, y){
    # circos.rect(xleft = x-1, ybottom = rep(0, length(y)), xright = x, ytop = 5)
    circos.points(x, 0.5)
  },
  track.height = 0.05
)

circos.clear()

dev.off()

# plot()
################






