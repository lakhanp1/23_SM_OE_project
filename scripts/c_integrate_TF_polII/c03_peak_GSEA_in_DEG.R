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

analysisName <- "peakset_enrichmet_in_DEG"
outDir <- here::here("analysis", "10_TF_polII_integration", "peak_GSEA_in_DEGs")
outPrefix <- paste(outDir, "/", analysisName, sep = "")
gseaPtOutDir <- paste(outDir, "/GSEA_plots", sep = "")

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

if(!dir.exists(gseaPtOutDir)){
  dir.create(path = gseaPtOutDir)
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

ptGseaLineDf <- NULL
ptGseaGeneDf <- NULL
ptGseaYlimDf <- NULL

rowId <- 1

for (rowId in 1:nrow(productionData)) {
  
  smTf <- productionData$SMTF[rowId]
  tfId <- productionData$tfId[rowId]
  degId <- productionData$degId[rowId]
  smTfName <- productionData$SMTF_name[rowId]
  
  ## prepare ranked polII DEG list
  degs <- dplyr::filter(combinedDegs, comparison == degId) %>% 
    dplyr::mutate(
      fpkmFilterSign = dplyr::if_else(
        condition = maxFpkm >= cutoff_fpkm, true = 1 * sign(shrinkLog2FC),
        false = 0, missing = 0
      ),
      rankMetric = -log10(pvalue) * sign(shrinkLog2FC)
    ) %>%
    dplyr::arrange(desc(rankMetric)) %>%
    dplyr::filter(!is.na(rankMetric))
  # dplyr::mutate(rankMetric = (-log10(pvalue) * sign(shrinkLog2FC))) %>%
  # dplyr::arrange(desc(rankMetric)) %>%
  # dplyr::filter(!is.na(rankMetric))
  
  geneList <- dplyr::select(degs, geneId, rankMetric) %>% 
    tibble::deframe()
  
  ## replace +Inf and -Inf values with max and min
  geneList[is.infinite(geneList) & geneList > 0] <- max(geneList[is.finite(geneList)]) + 1
  
  geneList[is.infinite(geneList) & geneList < 0] <- min(geneList[is.finite(geneList)]) - 1
  
  # barplot(geneList)
  
  ## prepare TF target gene list
  peakAn <- suppressMessages(readr::read_tsv(file = tfInfoList[[tfId]]$peakAnno)) %>% 
    dplyr::filter(peakPval >= cutoff_macs2Pval) %>% 
    dplyr::filter(!peakCategory %in% c("intergenic", "blacklist")) %>% 
    dplyr::filter(peakChr != "mito_A_nidulans_FGSC_A4")
  
  # table(peakAn$peakCategory)
  peakTargets <- list(unique(peakAn$geneId)) %>% purrr::set_names(nm = smTf)
  
  if(length(peakTargets[[smTf]]) == 0){
    next
  }
  
  cat("Running GSEA for ", smTf,"\n")
  
  fgseaRes <- fgsea::fgsea(
    pathways = peakTargets,
    stats = geneList,
    eps = 0,
    nPermSimple = 10000
  )
  
  fgseaRes$SM_TF <- smTf
  fgseaRes$SMTF_name <- smTfName
  
  peakGseaRes <- data.table::rbindlist(list(peakGseaRes, fgseaRes))
  
  ptData <- gsea_plot_data(geneset = peakTargets[[smTf]], stats = geneList)
  
  ptGseaLineDf <- dplyr::mutate(
    ptData$df,
    SM_TF = smTf,
    SMTF_name = smTfName,
    significance = if_else(
      condition = fgseaRes$pval <= 0.05,
      true = "significant", false = "non-significant"
    )
  ) %>% 
    dplyr::bind_rows(ptGseaLineDf)
  
  ptGseaGeneDf <- dplyr::mutate(
    ptData$genes, SM_TF = smTf, SMTF_name = smTfName,
  ) %>% 
    dplyr::bind_rows(ptGseaGeneDf)
  
  ptGseaYlimDf <- tibble::tibble(
    top = ptData$top, bottom = ptData$bottom,
    SM_TF = smTf, SMTF_name = smTfName,
  ) %>% 
    dplyr::bind_rows(ptGseaYlimDf)
  
  
  # plotGseaTable(pathways = peakTargets, stats = geneList, fgseaRes = fgseaRes)
  
  pt_gsea <- fgsea::plotEnrichment(
    pathway = peakTargets[[smTf]],
    stats = geneList,
    ticksSize = 0.1
  ) +
    labs(
      title = paste(smTfName, "peak targets' GSEA on OE/WT DEG list"),
      subtitle = paste(
        "NES =", round(fgseaRes$NES, digits = 4),
        "|| pval =", format(fgseaRes$pval, digits = 5)
      ),
      x = "Rank",
      y = "Enrichment Score"
    ) +
    ptTheme
  
  gseaPtOutPrefix <- paste(gseaPtOutDir, "/", smTf, ".GSEA_plot", sep = "")
  
  png(filename = paste(gseaPtOutPrefix, ".fgsea_plot.png", sep = ""),
      width = 3500, height = 1500, res = 400)

  print(pt_gsea)

  dev.off()
  
}

peakGseaResDf <- as.data.frame(peakGseaRes) %>% 
  dplyr::mutate(
    leadingEdgeLen = purrr::map_dbl(.x = leadingEdge, .f = length),
    leadingEdge = purrr::map_chr(.x = leadingEdge, .f = ~ paste(.x, collapse = ";"))
  ) %>% 
  dplyr::select(SM_TF, SMTF_name, -leadingEdge, everything(), leadingEdge)

readr::write_tsv(x = peakGseaResDf, file = paste(outPrefix, ".fgsea.tab", sep = ""))

peakGseaResDf <- dplyr::arrange(peakGseaResDf, NES) %>% 
  dplyr::mutate(
    SMTF_name = forcats::as_factor(SMTF_name),
    log10Padj = -log10(padj),
    sig = if_else(condition = padj <= 0.05, true = "**", false = "-", missing = "GSEA failed")
  )



pt_nes <- ggplot(data = peakGseaResDf, mapping = aes(x = NES, y = SMTF_name)) +
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

##################################################################################

ptGseaLineDf <- dplyr::mutate(
  ptGseaLineDf,
  SMTF_name = forcats::fct_relevel(.f = SMTF_name, levels(peakGseaResDf$SMTF_name))
)

ptGseaYlimDf <- dplyr::mutate(
  ptGseaYlimDf,
  SMTF_name = forcats::fct_relevel(.f = SMTF_name, levels(peakGseaResDf$SMTF_name))
)

ptGseaGeneDf <- dplyr::mutate(
  ptGseaGeneDf,
  SMTF_name = forcats::fct_relevel(.f = SMTF_name, levels(peakGseaResDf$SMTF_name))
)

pt_combined <- ggplot(data = ptGseaLineDf, mapping = aes(x = x, y = y)) +
  ggbeeswarm::geom_beeswarm(
    data = ptGseaGeneDf, mapping = aes(x = x, y = 0),
    groupOnX = FALSE
  ) +
  geom_point(mapping = aes(color = significance), size = 0.1) +
  geom_line(mapping = aes(color = significance)) +
  geom_hline(
    data = ptGseaYlimDf, mapping = aes(yintercept = top),
    colour = "red", linetype = "dashed"
  ) +
  geom_hline(
    data = ptGseaYlimDf, mapping = aes(yintercept = bottom),
    colour = "red", linetype = "dashed"
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_color_manual(
    values = c("significant" = "green", "non-significant" = "blue")
  ) +
  facet_wrap(
    facets = vars(SMTF_name), ncol = 4,
    scales = "free_y", strip.position = "right", dir = "h"
  ) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    legend.position = c(1, 0),
    legend.direction = "horizontal",
    legend.justification = c(1, 0),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank()
  )

png(filename = paste(outPrefix, ".fgsea_combined.png", sep = ""),
    width = 7000, height = 4000, res = 500)

print(pt_combined)

dev.off()


##################################################################################




