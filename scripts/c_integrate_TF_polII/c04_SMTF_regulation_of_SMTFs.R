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

outDir <- here::here("analysis", "10_TF_polII_integration", "OESMTF_self_regulation")
outPrefix <- paste(outDir, "/", sep = "")

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_polIIDegs <- here::here("analysis", "06_polII_diff", "polII_DEGs.fpkm_filtered.combined.tab")

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

genesDf <- as.data.frame(GenomicFeatures::genes(x = txDb)) %>% 
  dplyr::select(geneId = gene_id, strand)

##################################################################################

smTfs <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::select(SMTF) %>% 
  dplyr::distinct()

rowId <- 1

mergedData <- NULL

for (rowId in 1:nrow(productionData)) {
  
  tfSampleId <- productionData$tfId[rowId]
  degId <- productionData$degId[rowId]
  
  ## extract peak annotation
  peakAn <- suppressMessages(readr::read_tsv(tfInfoList[[tfSampleId]]$peakAnno)) %>% 
    dplyr::mutate(
      tfSampleId = tfSampleId,
    ) %>% 
    dplyr::filter(peakPval >= cutoff_macs2Pval) %>% 
    dplyr::select(peakId, tfSampleId, peakPval, geneId, peakAnnotation, peakDist, summitDist)
  
  
  ## extract DEG data
  diffData <-dplyr::filter(combinedDegs, comparison == degId) %>% 
    dplyr::select(
      geneId, comparison, !!col_lfc, shrinkLog2FC,
      pvalue, padj, maxFpkm, fpkmFilter
    ) %>% 
    dplyr::mutate(comparison = degId)
  
  bindingDegData <- dplyr::left_join(
    x = smTfs, y = diffData, by = "geneId"
  ) %>% 
    dplyr::left_join(y = peakAn, by = "geneId") %>% 
    dplyr::mutate(
      OESMTF = !!productionData$SMTF[rowId],
      OESMTF_name = !!productionData$SMTF_name[rowId]
    ) %>% 
    dplyr::group_by(geneId) %>% 
    dplyr::arrange(desc(peakPval), .by_group = TRUE) %>% 
    dplyr::slice(1L) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(
      OESMTF, geneId, everything()
    )
  
  mergedData <- dplyr::bind_rows(mergedData, bindingDegData)
  
}

mergedData2 <- dplyr::mutate(
  mergedData,
  significance = dplyr::case_when(
    !!sym(col_pval) <= cutoff_fdr & fpkmFilter == "pass" ~ "significant",
    TRUE ~ "non-significant"
  ),
  selfBinding = dplyr::case_when(
    !is.na(peakId) & OESMTF == geneId ~ "self",
    !is.na(peakId) & OESMTF != geneId ~ "cross",
    TRUE ~ "no-binding"
  ),
  selfBinding = forcats::fct_relevel(.f = selfBinding, "self")
) %>% 
  dplyr::arrange(selfBinding)

mergedData2$geneName <- AnnotationDbi::mapIds(
  x = orgDb, keys = mergedData2$geneId, column = "GENE_NAME", keytype = "GID"
)


## clusters the data to determine the row and column orders
dataMat <- mergedData2 %>% 
  dplyr::mutate(
    matrixVal = dplyr::case_when(
      significance == "significant" & !is.na(peakId) ~ 3*sign(log2FoldChange),
      significance == "significant" & is.na(peakId) ~ 2*sign(log2FoldChange),
      significance == "non-significant" & !is.na(peakId) ~ 1,
      TRUE ~ 0
    )
  ) %>% 
  tidyr::pivot_wider(
    id_cols = c(OESMTF_name),
    names_from = geneName,
    values_from = matrixVal
  ) %>% 
  tibble::column_to_rownames(var = "OESMTF_name") %>% 
  as.matrix()

hc_oeTf <- hclust(d = dist(x = t(dataMat[, rownames(dataMat)])))

oetf_levels <- hc_oeTf$labels[hc_oeTf$order]


hc_nonoeTf <- hclust(
  d = dist(t(dataMat[, setdiff(x = colnames(dataMat), y = rownames(dataMat))]))
)

smtf_levels <- c(oetf_levels, hc_nonoeTf$labels[hc_nonoeTf$order])
  
# smtf_levels <- intersect(mergedData2$geneName, mergedData2$OESMTF_name)
# oetf_levels <- intersect(mergedData2$geneName, mergedData2$OESMTF_name)

mergedData2 <- dplyr::mutate(
  mergedData2,
  geneName = forcats::fct_relevel(.f = geneName, smtf_levels),
  OESMTF_name = forcats::fct_relevel(.f = OESMTF_name, oetf_levels)
)

##################################################################################

pt_binding <- ggplot(data = mergedData2, mapping = aes(x = geneName, y = OESMTF_name)) +
  geom_point(
    mapping = aes(color = selfBinding), size = 3
  ) +
  scale_color_manual(
    values = c("self" = "green", "cross" = "black", "no-binding" = alpha("white", alpha = 0)),
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

ggdraw(aligned_plots[[2]])

ggsave(
  filename = paste(outPrefix, "OESMTF_binding.SMTFs.png", sep = ""), plot = ggdraw(aligned_plots[[1]]),
  width = 16, height = 10
)

ggsave(
  filename = paste(outPrefix, "OESMTF_binding_DEG.SMTFs.png", sep = ""), plot = ggdraw(aligned_plots[[2]]),
  width = 16, height = 10
)








