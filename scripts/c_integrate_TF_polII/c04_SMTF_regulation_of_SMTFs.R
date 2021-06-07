suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(gginnards))

## This script plots the SMTF's binding and log2(OE/WT) data for all the SMTFs as matrix

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

##################################################################################

outDir <- here::here("analysis", "10_TF_polII_integration", "OESMTF_binding_DEG")
outPrefix <- paste(outDir, "/", sep = "")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_dataSummary <- here::here("data", "reference_data", "production_data.summary.tab")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")

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

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

dataSummary <- suppressMessages(readr::read_tsv(file = file_dataSummary)) %>% 
  dplyr::filter(has_TF_ChIP == "has_data", has_polII_ChIP == "has_data", copyNumber == "sCopy") %>% 
  dplyr::arrange(SM_ID, geneId)


tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = dataSummary$tfId,
  dataPath = TF_dataPath)

tfInfoList <- purrr::transpose(tfInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% dataSummary$degId)

rnaseqInfoList <- purrr::transpose(rnaseqInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "comparison"))

genesDf <- as.data.frame(GenomicFeatures::genes(x = txDb)) %>% 
  dplyr::select(geneId = gene_id, strand)

##################################################################################

smTfs <- AnnotationDbi::select(
  x = orgDb, keys = keys(x = orgDb, keytype = "SM_CLUSTER"),
  columns = c("SM_ID", "GID", "TF_GENE"), keytype = "SM_CLUSTER"
) %>% 
  dplyr::rename(geneId = GID) %>% 
  dplyr::filter(!is.na(TF_GENE))

rowId <- 1

mergedData <- NULL

for (rowId in 1:nrow(dataSummary)) {
  
  tfSampleId <- dataSummary$tfId[rowId]
  degId <- dataSummary$degId[rowId]
  
  ## extract peak annotation
  peakAn <- suppressMessages(readr::read_tsv(tfInfoList[[tfSampleId]]$peakAnno)) %>% 
    dplyr::mutate(
      tfSampleId = tfSampleId,
    ) %>% 
    dplyr::filter(peakPval >= cutoff_macs2Pval)
  
  
  ## extract DEG data
  diffData <- suppressMessages(readr::read_tsv(file = rnaseqInfoList[[degId]]$deg)) %>% 
    dplyr::select(geneId, !!col_lfc, shrinkLog2FC, pvalue, padj) %>% 
    dplyr::mutate(comparison = degId)
  
  bindingDegData <- dplyr::left_join(
    x = diffData, y = peakAn, by = "geneId"
  ) %>% 
    dplyr::mutate(
      OESMTF = !!tfInfoList[[tfSampleId]]$SM_TF
    ) %>% 
    dplyr::left_join(
      y = smTfs, by = "geneId"
    ) %>% 
    dplyr::filter(!is.na(SM_ID))
  
  mergedData <- dplyr::bind_rows(mergedData, bindingDegData)
  
}

mergedData2 <- dplyr::mutate(
  mergedData,
  significance = if_else(
    condition = !!sym(col_pval) <= cutoff_fdr, true = "significant",
    false = "non-significant", missing = "non-significant"
  ),
  selfBinding = dplyr::if_else(
    condition = OESMTF == geneId, true = "self", false = "cross"
  ),
  selfBinding = dplyr::if_else(
    condition = is.na(peakId), true = "no-binding", false = selfBinding
  ),
  selfBinding = forcats::fct_relevel(.f = selfBinding, "self")
) %>% 
  dplyr::arrange(selfBinding)

mergedData2$geneName <- AnnotationDbi::mapIds(
  x = orgDb, keys = mergedData2$geneId, column = "GENE_NAME", keytype = "GID"
)

mergedData2$OESMTF_name <- AnnotationDbi::mapIds(
  x = orgDb, keys = mergedData2$OESMTF, column = "GENE_NAME", keytype = "GID"
)


fctLevels <- intersect(mergedData2$geneName, mergedData2$OESMTF_name)

mergedData2 <- dplyr::mutate(
  mergedData2,
  geneName = forcats::fct_relevel(.f = geneName, fctLevels),
  OESMTF_name = forcats::fct_relevel(.f = OESMTF_name, fctLevels)
)

pt_binding <- ggplot(data = mergedData2, mapping = aes(x = geneName, y = OESMTF_name)) +
  geom_point(
    mapping = aes(color = selfBinding), size = 3
  ) +
  scale_color_manual(
    values = c("self" = "green", "cross" = "black", "no-binding" = NA),
    breaks = c("self", "cross"), name = "TF binding"
    
  ) +
  labs(
    title = "Self and cross binding of SMTFs based on TF ChIPseq",
    x = "TF genes",
    y = "SMTF-OE strain"
  ) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 13, face = "bold"),
    panel.grid = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 15, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold")
  )

pt_bindingLfc <- pt_binding +
  geom_tile(
    mapping = aes(fill = log2FoldChange, alpha = significance)
  ) +
  scale_fill_gradient2(
    name = paste("RNA-polII log2(", "OE/WT", ")", sep = ""),
    low = "#313695", mid = "#F7F7F7", high = "#d73027", midpoint = 0,
    limit = c(-2.5, 2.5), oob = scales::squish
  ) +
  scale_alpha_manual(
    values = c("significant" = 1, "non-significant" = 0),
    breaks = NULL
  ) +
  labs(
    title = "Self and cross regulation of SMTFs based on TF binding and transcription response",
    x = "TF genes",
    y = "SMTF-OE strain"
  )

pt_bindingLfc <- gginnards::move_layers(x = pt_bindingLfc, match_type = "GeomPoint", position = "top")

aligned_plots <- cowplot::align_plots(plotlist = list(pt_binding, pt_bindingLfc), align = "v")

ggsave(
  filename = paste(outPrefix, "OESMTF_binding.SMTFs.png", sep = ""), plot = ggdraw(aligned_plots[[1]]),
  width = 16, height = 10
)

ggsave(
  filename = paste(outPrefix, "OESMTF_binding_DEG.SMTFs.png", sep = ""), plot = ggdraw(aligned_plots[[2]]),
  width = 16, height = 10
)








