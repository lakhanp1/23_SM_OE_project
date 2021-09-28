suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(gginnards))

## this script plots SMTFOE strain's respective SM cluster genes and shows:
## 1) log2(OE/WT) 
## 2) binding status by its respective SMTF

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

##################################################################################

analysisName <- "OESMTF_binding_DEG.self_clusters"
outDir <- here::here("analysis", "10_TF_polII_integration", "OESMTF_self_regulation")
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

cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1

col_lfc <- "log2FoldChange"
col_pval <- "pvalue"
##################################################################################

productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::filter(has_polII_ChIP == "has_data", has_TF_ChIP == "has_data", copyNumber == "sCopy") %>% 
  dplyr::filter(!is.na(SM_ID)) %>% 
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

genesDf <- as.data.frame(GenomicFeatures::genes(x = txDb)) %>% 
  dplyr::select(geneId = gene_id, strand)

genePos <- as.data.frame(GenomicFeatures::genes(x = txDb)) %>% 
  dplyr::select(geneId = gene_id, chr = seqnames, start, end, strand)

combinedDegs <- suppressMessages(readr::read_tsv(file = file_polIIDegs))


##################################################################################

smGenes <- AnnotationDbi::select(
  x = orgDb, keys = keys(x = orgDb, keytype = "SM_CLUSTER"),
  columns = c("SM_ID", "GID", "TF_GENE"), keytype = "SM_CLUSTER"
) %>% 
  dplyr::rename(geneId = GID) %>% 
  dplyr::left_join(y = genePos, by = "geneId") %>% 
  dplyr::group_by(SM_ID) %>% 
  dplyr::arrange(chr, start, .by_group = TRUE) %>% 
  dplyr::mutate(index = row_number()) %>% 
  dplyr::ungroup()

rowId <- 1

mergedData <- NULL

for (rowId in 1:nrow(productionData)) {
  
  tfSampleId <- productionData$tfId[rowId]
  degId <- productionData$degId[rowId]
  
  clusterGenes <- dplyr::filter(smGenes, SM_ID == productionData$SM_ID[rowId])
  
  ## extract peak annotation
  peakAn <- suppressMessages(readr::read_tsv(tfInfoList[[tfSampleId]]$peakAnno)) %>% 
    dplyr::filter(peakPval >= cutoff_macs2Pval)
  
  ## extract DEG data
  diffData <- dplyr::filter(combinedDegs, comparison == degId) %>% 
    dplyr::select(geneId, !!col_lfc, shrinkLog2FC, pvalue, padj, maxFpkm, fpkmFilter)
  
  bindingDegData <- dplyr::left_join(
    x = clusterGenes, y = diffData, by = "geneId"
  ) %>% 
    dplyr::left_join(y = peakAn, by = "geneId") %>% 
    dplyr::mutate(
      OESMTF = !!productionData$SMTF[rowId],
      OESMTF_name = !!productionData$SMTF_name[rowId]
    )
  
  mergedData <- dplyr::bind_rows(mergedData, bindingDegData)
  
}


mergedData2 <- dplyr::mutate(
  mergedData,
  significance = dplyr::case_when(
    !!sym(col_pval) <= cutoff_fdr & fpkmFilter == "pass" ~ "significant",
    TRUE ~ "non-significant"
  ),
  binding = if_else(
    condition = !is.na(peakId), true = "bound", false = "not-bound"
  ),
  tfGene = if_else(
    condition = !is.na(TF_GENE), true = "TF", false = "non-TF"
  )
) %>% 
  tidyr::unite(col = "tf_cluster_grp", OESMTF_name, SM_CLUSTER, sep = ": ")


pltTitle <- "RNA-polII ChIPseq log2(SMTF-OE/WT) of SM clusters in its own TF overexpression data and SMTF ChIPseq binding"

pt_bindingLfc <- ggplot(
  data = mergedData2,
) +
  geom_tile(
    mapping = aes(x = index, y = tf_cluster_grp, fill = log2FoldChange, alpha = significance),
    size = 0.2, height = 1, color = "black") +
  geom_point(
    mapping = aes(x = index , y = "peak", color = binding, shape = tfGene),
    size = 3
  ) +
  scale_fill_gradient2(
    name = paste("polII-ChIPseq log2(", "OE/WT", ")", sep = ""),
    low = "#313695", mid = "#F7F7F7", high = "#d73027", midpoint = 0,
    limit = c(-2.5, 2.5), oob = scales::squish
  ) +
  scale_alpha_manual(
    values = c("significant" = 1, "non-significant" = 0),
    breaks = NULL
  ) +
  scale_colour_manual(
    name = "Peak", breaks = "bound",
    values = c("bound" = "black", "not-bound" = alpha("white", alpha = 0))
  ) +
  scale_shape_manual(
    name = "Gene type", breaks = c("TF", "non-TF"),
    values = c("TF" = 17, "non-TF" = 16)
  ) +
  scale_x_continuous(expand = expansion(add = c(0.0, 0.0))) +
  scale_y_discrete(
    expand = expansion(add = c(0.0, 0.5))
  ) +
  facet_wrap(
    facets = . ~ tf_cluster_grp, scales = "free_y",
    ncol = 4, dir = "v"
  ) +
  ggtitle(pltTitle) +
  guides(
    color = guide_legend(override.aes = list(shape = 17, size = 5)),
    shape = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.15, "lines"),
    strip.background = element_rect(fill="white", size = 0.2),
    strip.text.x = element_text(size = 12, hjust = 0),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    plot.margin = unit(rep(0.2, 4), "cm")
  )


ggsave(filename = paste(outPrefix, ".tiles.png", sep = ""),
       plot = pt_bindingLfc, width = 14, height = 8)







