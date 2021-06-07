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
##################################################################################

analysisName <- "peakset_enrichmet_in_DEG"
outDir <- here::here("analysis", "10_TF_polII_integration", "peakset_enrichmet_in_DEG")
outPrefix <- paste(outDir, "/", analysisName, sep = "")
gseaPtOutDir <- paste(outDir, "/GSEA_plots", sep = "")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_dataSummary <- here::here("data", "reference_data", "raw_data_summary.tab")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")

file_peakSummary <- here::here("analysis", "02_QC_TF", "TF_ChIP_summary.best_replicates.tab")
file_degSummary <- here::here("analysis", "08_polII_diff_downstream", "01_polII_DEGs_summary", "polII_DEGs_summary.stats.tab")

TF_dataPath <- here::here("data", "TF_data")
diffDataPath <- here::here("analysis", "06_polII_diff")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

cutoff_macs2Pval <- 20

col_lfc <- "log2FoldChange"
cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1

##################################################################################

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

if(!dir.exists(gseaPtOutDir)){
  dir.create(path = gseaPtOutDir)
}

dataSummary <- suppressMessages(readr::read_tsv(file = file_dataSummary)) %>% 
  dplyr::filter(has_TF_ChIP == "has_data", has_polII_ChIP == "has_data")

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath)

rnaseqInfoList <- purrr::transpose(rnaseqInfo) %>% 
  purrr::set_names(nm = purrr::map(., "comparison"))

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = dataSummary$tfId,
  dataPath = TF_dataPath)

tfInfoList <- purrr::transpose(tfInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))


##################################################################################


ptTheme <- theme_bw() +
  theme(
    axis.text = element_text(size = 16, face = "bold"),
    title = element_text(size = 16, face = "bold"),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    plot.margin = margin(0.25, 0.75, 0.25, 0.25, "cm"),
    panel.grid = element_blank()
  )

masterCombDf <- NULL
peakGseaRes <- NULL

rowId <- 21

for (rowId in 1:nrow(dataSummary)) {
  
  smTf <- dataSummary$geneId[rowId]
  tfId <- dataSummary$tfId[rowId]
  degId <- dataSummary$degId[rowId]
  
  ## prepare ranked polII DEG list
  degs <- suppressMessages(readr::read_tsv(file = rnaseqInfoList[[degId]]$deg)) %>% 
    dplyr::mutate(rankMetric = (-log10(pvalue) * sign(shrinkLog2FC))) %>%
    dplyr::arrange(desc(rankMetric)) %>% 
    dplyr::filter(!is.na(rankMetric))
  
  downDegs <- dplyr::filter(degs, padj <= cutoff_fdr & !!sym(col_lfc) <= cutoff_down) %>% 
    dplyr::mutate(category = "down")
  
  upDegs <- dplyr::filter(degs, padj <= cutoff_fdr & !!sym(col_lfc) >= cutoff_up) %>% 
    dplyr::mutate(category = "up")
  
  geneList <- dplyr::select(degs, geneId, rankMetric) %>% 
    tibble::deframe()
  
  ## replace +Inf and -Inf values with max and min
  geneList[is.infinite(geneList) & geneList > 0] <- max(geneList[is.finite(geneList)]) + 1
  
  geneList[is.infinite(geneList) & geneList < 0] <- min(geneList[is.finite(geneList)]) - 1
  
  ## prepare TF target gene list
  peakAn <- suppressMessages(readr::read_tsv(file = tfInfoList[[tfId]]$peakAnno)) %>% 
    dplyr::filter(peakPval >= cutoff_macs2Pval) %>% 
    dplyr::filter(!peakCategory %in% c("intergenic", "blacklist"))
  
  # table(peakAn$peakCategory)
  peakTargets <- unique(peakAn$geneId)
  
  cat("Running GSEA for ", smTf,"\n")
  
  # fgseaRes <- fgsea::fgsea(
  #   pathways = list(peakTargets) %>% purrr::set_names(nm = smTf),
  #   stats = geneList,
  #   eps = 0.0
  # )
  # 
  # peakGseaRes <- data.table::rbindlist(list(peakGseaRes, fgseaRes))
  
  pt_gsea <- plotEnrichment(
    pathway = peakTargets,
    stats = geneList,
    ticksSize = 0.1
  ) +
    labs(
      title = paste(smTf, "bound gene set GSEA on ranked DEG list"),
      x = "Rank",
      y = "Enrichment Score"
    ) +
    ptTheme
  
  gseaPtOutPrefix <- paste(gseaPtOutDir, "/", smTf, ".GSEA_plot", sep = "")
  
  png(filename = paste(gseaPtOutPrefix, ".fgsea_plot.png", sep = ""),
      width = 3000, height = 1500, res = 400)
  print(pt_gsea)
  dev.off()
  
  
}


peakGseaRes <- as.data.frame(peakGseaRes) %>% 
  dplyr::mutate(
    leadingEdgeLen = purrr::map_dbl(.x = leadingEdge, .f = length),
    leadingEdge = purrr::map_chr(.x = leadingEdge, .f = ~ paste(.x, collapse = ";"))
  ) %>% 
  dplyr::select(-leadingEdge, everything(), leadingEdge)

readr::write_tsv(x = peakGseaRes, path = paste(outPrefix, ".fgsea.tab", sep = ""))

peakGseaRes <- dplyr::arrange(peakGseaRes, NES) %>% 
  dplyr::mutate(
    pathway = forcats::as_factor(pathway),
    log10Padj = -log10(padj),
    sig = if_else(condition = padj <= 0.05, true = "**", false = "-", missing = "GSEA failed")
  )




pt_nes <- ggplot(data = peakGseaRes, mapping = aes(x = NES, y = pathway)) +
  geom_bar(mapping = aes(fill = sig), stat = "identity") +
  scale_fill_manual(
    name = "Significance",
    values = c("**" = "green", "-" = "black")
  ) +
  geom_text(mapping = aes(x = (0.2*sign(NES)*-1), label = round(log10Padj, 2))) +
  labs(
    title = "SM TFOE peak targets GSEA in polII DEG list: NES scores",
    x = "Normalized Enrichment Score"
  ) +
  ptTheme +
  theme(
    legend.position = c(0.8, 0.1),
    axis.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold"),
  )

png(filename = paste(outPrefix, ".fgsea_NES_bar.png", sep = ""), width = 2500, height = 3000, res = 300)
pt_nes
dev.off()






