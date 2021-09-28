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
## work on this to use new stats

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")
##################################################################################

analysisName <- "peak_DEG_stats"
outDir <- here::here("analysis", "10_TF_polII_integration", "stats")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_polIIDegs <- here::here("analysis", "06_polII_diff", "polII_DEGs.fpkm_filtered.combined.tab")

TF_dataPath <- here::here("data", "TF_data")
diffDataPath <- here::here("analysis", "06_polII_diff")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

cutoff_macs2Pval <- 20

cutoff_fpkm <- 10
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

productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::filter(has_polII_ChIP == "has_data", has_TF_ChIP == "has_data", copyNumber == "sCopy") %>% 
  dplyr::arrange(SM_ID, SMTF)

productionDataList <- purrr::transpose(productionData) %>% 
  purrr::set_names(nm = purrr::map(., "SMTF"))

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


##################################################################################

# masterSet <- list(
#   up = character(), down = character(), peaks = character()
# )

masterCombDf <- NULL

rowId <- 1

for (rowId in 1:nrow(productionData)) {
  
  smTf <- productionData$SMTF[rowId]
  tfId <- productionData$tfId[rowId]
  degId <- productionData$degId[rowId]
  smTfName <- productionData$SMTF_name[rowId]
  
  ## prepare ranked polII DEG list
  degs <- dplyr::filter(combinedDegs, comparison == degId) %>%
    dplyr::mutate(rankMetric = (-log10(pvalue) * sign(shrinkLog2FC))) %>%
    dplyr::arrange(desc(rankMetric)) %>% 
    dplyr::filter(!is.na(rankMetric))
  
  downDegs <- dplyr::filter(
    degs, padj <= cutoff_fdr & !!sym(col_lfc) <= cutoff_down, maxFpkm >= cutoff_fpkm
  ) %>% 
    dplyr::mutate(category = "down")
  
  upDegs <- dplyr::filter(
    degs, padj <= cutoff_fdr & !!sym(col_lfc) >= cutoff_up, maxFpkm >= cutoff_fpkm
  ) %>% 
    dplyr::mutate(category = "up")
  
  ## extract peak annotation
  peakAn <- suppressMessages(readr::read_tsv(tfInfoList[[tfId]]$peakAnno)) %>%
    dplyr::filter(
      peakPval >= cutoff_macs2Pval,
      !peakCategory %in% c("intergenic", "blacklist")
    ) %>%
    dplyr::mutate(hasPeak = TRUE) %>%
    dplyr::group_by(geneId) %>%
    dplyr::arrange(desc(peakPval), .by_group = TRUE) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::select(geneId, peakPval, peakAnnotation, peakCategory, peakPosition, hasPeak)
  
  # table(peakAn$peakCategory)
  peakTargets <- unique(peakAn$geneId)
  
  combDf <- tibble::tibble(
    SMTF = smTf,
    SMTF_name = smTfName,
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
  #   SMTF = smTf,
  #   comb = names(cmSize),
  #   size = cmSize,
  #   sizeFraction = round(x = 100*cmSize/sum(cmSize), digits = 4)
  # )
  
  masterCombDf <- dplyr::bind_rows(masterCombDf, combDf)
  
  
}

masterCombDf <- dplyr::arrange(masterCombDf, desc(peaks)) %>% 
  dplyr::mutate(
    SMTF = forcats::as_factor(SMTF),
    SMTF_name = forcats::as_factor(SMTF_name)
  )


degPeakOvlpDf <- dplyr::select(
  masterCombDf, SMTF, SMTF_name, up_peak, up_noPeak, down_peak, down_noPeak
) %>% 
  tidyr::pivot_longer(
    cols = !c(SMTF, SMTF_name),
    names_to = c("DEG", "with_peak"),
    names_pattern = "(\\w+)_(\\w+)",
    values_to = "count"
  ) %>% 
  dplyr::mutate(
    DEG = forcats::fct_relevel(DEG, "up"),
    with_peak = forcats::fct_relevel(with_peak, "peak")
  )

peakDegOvlpDf <- dplyr::select(
  masterCombDf, SMTF, SMTF_name, peaks, up_peak, down_peak, noDeg_peak
) %>% 
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
    plot.title = element_text(size = 14, face = "bold", hjust = 1),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    plot.margin = margin(0.25, 0.75, 0.25, 0.25, "cm"),
    panel.grid = element_blank()
  )

## peak count bar plot
pt_peakCount <- ggplot(data = masterCombDf) +
  geom_bar(mapping = aes(y = forcats::fct_rev(SMTF_name), x = peaks),
           stat = "identity", fill = col_peak, color = "black", width = 0.8) +
  scale_x_continuous(expand = expansion(add = 0)) +
  labs(title = "Peak count") +
  ptTheme

## DEG count bar plot
pt_deg <- dplyr::mutate(masterCombDf, down = -1*down) %>% 
  ggplot(mapping = aes(y = forcats::fct_rev(SMTF_name))) +
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
    mapping = aes(y = forcats::fct_rev(SMTF_name), x = count, fill = deg),
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
    mapping = aes(y = forcats::fct_rev(SMTF_name), x = count, fill = with_peak),
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
    mapping = aes(y = forcats::fct_rev(SMTF_name), x = count, fill = with_peak),
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


pt_mergedAll <- ggarrange(
  pt_peakCount, pt_peakOvlp, pt_downOvlp, pt_deg, pt_upOvlp,
  widths = c(1, 1), ncol = 5, align = "h",
  common.legend = TRUE, legend = "right"
)

# png(filename = paste(outPrefix, ".overlap.png", sep = ""), width = 6000, height = 3500, res = 320)
# # pdf(file = paste(outPrefix, ".overlap.pdf", sep = ""), width = 20, height = 10)
# pt_mergedAll
# dev.off()

ggsave(
  filename = paste(outPrefix, ".overlap.png", sep = ""),
  plot = pt_mergedAll, width = 18, height = 10
)


pt_peakDeg <- ggarrange(
  pt_peakCount, pt_deg,
  widths = c(1, 1), ncol = 2, align = "h",
  common.legend = TRUE, legend = "right"
)


ggsave(filename = paste(outPrefix, ".png", sep = ""), plot = pt_peakDeg, width = 12, height = 8)













