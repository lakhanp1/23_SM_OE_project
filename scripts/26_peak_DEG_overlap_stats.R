suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(fgsea))


## plot to show the peak and DEG overlap data 

rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")
##################################################################################

analysisName <- "peak_DEG_overlap"
outDir <- here::here("analysis", "10_TF_polII_integration", "peak_DEG_overlap")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_dataSummary <- here::here("data", "reference_data", "raw_data_summary.tab")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")

file_peakSummary <- here::here("analysis", "02_QC_TF", "TF_ChIP_summary.best_replicates.tab")
file_degSummary <- here::here("analysis", "06_polII_diff", "polII_DEG.stats.tab")

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

# masterSet <- list(
#   up = character(), down = character(), peaks = character()
# )

masterCombDf <- NULL

rowId <- 1

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
  
  combDf <- tibble::tibble(
    tf = smTf,
    up = length(upDegs$geneId),
    down = length(downDegs$geneId),
    nDegs = up + down,
    peaks = length(peakTargets),
    up_peak = length(intersect(upDegs$geneId, peakTargets)),
    up_noPeak = length(setdiff(upDegs$geneId, peakTargets)),
    down_peak = length(intersect(downDegs$geneId, peakTargets)),
    down_noPeak = length(setdiff(downDegs$geneId, peakTargets))
  ) %>% 
    dplyr::mutate(
      noDeg_peak = peaks - (up_peak + down_peak)
    )
  
  # ## Upset plot
  # geneSets <- list(
  #   up = upDegs$geneId, down = downDegs$geneId, peaks = peakTargets
  # )
  # 
  # masterSet$up <- union(masterSet$up, geneSets$up)
  # masterSet$down <- union(masterSet$down, geneSets$down)
  # masterSet$peaks <- union(masterSet$peaks, geneSets$peaks)
  # 
  # cm = make_comb_mat(geneSets)
  # # set_size(cm)
  # # comb_name(cm)
  # # comb_size(cm)
  # # comb_degree(cm)
  # # set_name(cm)
  # # t.comb_mat(cm)
  # cmSize <- comb_size(cm)
  # 
  # combDf <- tibble::tibble(
  #   tf = smTf,
  #   comb = names(cmSize),
  #   size = cmSize,
  #   sizeFraction = round(x = 100*cmSize/sum(cmSize), digits = 4)
  # )
  
  masterCombDf <- dplyr::bind_rows(masterCombDf, combDf)
  
  
}

masterCombDf <- dplyr::arrange(masterCombDf, desc(peaks)) %>% 
  dplyr::mutate(tf = forcats::as_factor(tf))


degPeakOvlpDf <- dplyr::select(masterCombDf, tf, up_peak, up_noPeak, down_peak, down_noPeak) %>% 
  tidyr::pivot_longer(
    cols = -tf,
    names_to = c("DEG", "with_peak"),
    names_pattern = "(\\w+)_(\\w+)",
    values_to = "count"
  ) %>% 
  dplyr::mutate(
    DEG = forcats::fct_relevel(DEG, "up"),
    with_peak = forcats::fct_relevel(with_peak, "peak")
  )

peakDegOvlpDf <- dplyr::select(masterCombDf, tf, peaks, up_peak, down_peak, noDeg_peak) %>% 
  tidyr::pivot_longer(
    cols = ends_with("_peak"),
    names_to = "deg",
    names_pattern = "(\\w+)_peak",
    values_to = "count"
  ) %>% 
  dplyr::mutate(
    deg = forcats::fct_relevel(deg, "down", "noDeg", "up")
  )

col_peak <- "#FFC20A"
col_down <- "#313695"
col_up <- "#a50026"

ptTheme <- theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.text.x = element_text(face = "bold"),
    axis.title = element_blank(),
    title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    plot.margin = margin(0.25, 0.75, 0.25, 0.25, "cm"),
    panel.grid = element_blank()
  )

## peak count bar plot
pt_peakCount <- ggplot(data = masterCombDf) +
  geom_bar(mapping = aes(y = forcats::fct_rev(tf), x = peaks),
           stat = "identity", fill = col_peak, color = "black", width = 0.8) +
  scale_x_continuous(expand = expansion(add = 0)) +
  labs(title = "Peak count") +
  ptTheme

## DEG count bar plot
pt_deg <- dplyr::mutate(masterCombDf, down = -1*down) %>% 
  ggplot(mapping = aes(y = forcats::fct_rev(tf))) +
  geom_bar(
    mapping = aes(x = up), stat = "identity", fill = col_up,
    width = 0.8, color = "black"
  ) +
  geom_bar(
    mapping = aes(x = down), stat = "identity", fill = col_down,
    width = 0.8, color = "black"
  ) +
  labs(title = "PolII DEG count") +
  ptTheme


## stacked bar plot: peak overlap with up and down DEGs
pt_peakOvlp <- ggplot(data = peakDegOvlpDf) +
  geom_bar(
    mapping = aes(y = forcats::fct_rev(tf), x = count, fill = deg),
    stat = "identity", position = position_fill(reverse = TRUE),
    width = 0.8, color = "black"
  ) +
  scale_fill_manual(
    breaks = NULL,
    values = c("down" = col_down, "noDeg" = col_peak, "up" = col_up)
  ) +
  scale_x_continuous(labels = scales::percent, expand = expansion(add = 0)) +
  labs(title = "% peaks targets as DEG") +
  ptTheme


## stacked bar plot: up DEG overlap with peak
pt_upOvlp <- dplyr::filter(degPeakOvlpDf, DEG == "up") %>% 
  ggplot() +
  geom_bar(
    mapping = aes(y = forcats::fct_rev(tf), x = count, fill = with_peak),
    stat = "identity", position = position_fill(reverse = TRUE),
    width = 0.8, color = "black"
  ) +
  scale_fill_manual(
    breaks = NULL,
    values = c("peak" = col_peak, "noPeak" = col_up)
  ) +
  scale_x_continuous(labels = scales::percent, expand = expansion(add = 0)) +
  labs(title = "% up DEG with peaks") +
  ptTheme

## stacked bar plot: down DEG overlp with peak
pt_downOvlp <- dplyr::filter(degPeakOvlpDf, DEG == "down") %>% 
  ggplot() +
  geom_bar(
    mapping = aes(y = forcats::fct_rev(tf), x = count, fill = with_peak),
    stat = "identity", position = position_fill(reverse = TRUE),
    width = 0.8, color = "black"
  ) +
  scale_fill_manual(
    breaks = NULL,
    values = c("peak" = col_peak, "noPeak" = col_down)
  ) +
  scale_x_continuous(labels = scales::percent, expand = expansion(add = 0)) +
  labs(title = "% down DEG with peaks") +
  ptTheme


pt_merged <- ggarrange(
  pt_peakCount, pt_peakOvlp, pt_downOvlp, pt_deg, pt_upOvlp,
  widths = c(1, 1), ncol = 5, align = "h",
  common.legend = TRUE, legend = "right"
)

png(filename = paste(outPrefix, ".stats.png", sep = ""), width = 6000, height = 3500, res = 320)
# pdf(file = paste(outPrefix, ".stats.pdf", sep = ""), width = 20, height = 10)
pt_merged
dev.off()




















