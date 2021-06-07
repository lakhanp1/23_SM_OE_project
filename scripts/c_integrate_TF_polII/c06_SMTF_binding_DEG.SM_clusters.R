suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(gginnards))

## this script plots tiles of SM clusters showing following data for each SMTF:
## 1) log2(OE/WT) 
## 2) binding status by its respective SMTF

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

##################################################################################

analysisName <- "OESMTF_binding_DEG.SM_clusters"
outDir <- here::here("analysis", "10_TF_polII_integration", "OESMTF_binding_DEG")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

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

dataSummary <- suppressMessages(readr::read_tsv(file = file_dataSummary)) %>% 
  dplyr::filter(has_TF_ChIP == "has_data", has_polII_ChIP == "has_data", copyNumber == "sCopy") %>% 
  dplyr::arrange(SM_ID, geneId)

dataSummary$OESMTF_name <- AnnotationDbi::mapIds(
  x = orgDb, keys = dataSummary$geneId, column = "GENE_NAME", keytype = "GID"
)

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

## extract all SM cluster genes
geneset <- AnnotationDbi::select(
  x = orgDb, keys = keys(x = orgDb, keytype = "SM_ID"),
  columns = c("GID", "GENE_NAME", "SM_CLUSTER", "TF_GENE"), keytype = "SM_ID"
) %>% 
  dplyr::rename(geneId = GID) %>% 
  dplyr::mutate(
    geneLabel = paste(geneId, "(", GENE_NAME, ")", sep = ""),
    geneLabel = if_else(condition = geneId == GENE_NAME, true = geneId, false = geneLabel)
  )

genePos <- as.data.frame(genes(txDb)) %>% 
  dplyr::select(geneId = gene_id, chr = seqnames, start, end, strand)

geneset <- dplyr::left_join(x = geneset, y = genePos, by = "geneId")

##################################################################################
## get DEG data for each SM_TF cluster from its own polII_DEG set
rowId <- 1
mergedData <- NULL


for (rowId in 1:nrow(dataSummary)) {
  
  tfSampleId <- dataSummary$tfId[rowId]
  degId <- dataSummary$degId[rowId]
  
  ## extract peak annotation
  peakAn <- suppressMessages(readr::read_tsv(tfInfoList[[tfSampleId]]$peakAnno)) %>% 
    dplyr::filter(peakPval >= cutoff_macs2Pval) %>% 
    dplyr::mutate(hasPeak = TRUE) %>% 
    dplyr::group_by(geneId) %>% 
    dplyr::arrange(desc(peakPval), .by_group = TRUE) %>% 
    dplyr::slice(1L) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(geneId, peakPval, peakAnnotation, peakCategory, peakPosition, hasPeak)
  
  ## extract DEG data
  diffData <- suppressMessages(readr::read_tsv(file = rnaseqInfoList[[degId]]$deg)) %>% 
    dplyr::select(geneId, !!col_lfc, shrinkLog2FC, pvalue, padj) %>% 
    dplyr::mutate(comparison = degId)
  
  
  bindingDegData <- dplyr::left_join(
    x = geneset, y = diffData, by = "geneId"
  ) %>% 
    dplyr::left_join(y = peakAn, by = "geneId") %>% 
    dplyr::mutate(
      OESMTF = !!dataSummary$geneId[rowId],
      OESMTF_name = !!dataSummary$OESMTF_name[rowId]
    ) %>% 
    dplyr::select(
      OESMTF, geneId, GENE_NAME, everything()
    )
  
  mergedData <- dplyr::bind_rows(mergedData, bindingDegData)
  
}


mergedData2 <- dplyr::mutate(
  mergedData,
  significance = if_else(
    condition = !!sym(col_pval) <= cutoff_fdr, true = "significant",
    false = "non-significant", missing = "non-significant"
  ),
  binding = dplyr::if_else(
    condition = hasPeak, true = "bound", false = "not-bound", missing = "not-bound"
  ),
  tfGene = if_else(
    condition = !is.na(TF_GENE), true = "TF", false = "non-TF"
  )
)

###########################################################################

pdf(file = paste(outPrefix, ".tiles.pdf", sep = ""), width = 15, height = 8, onefile = TRUE)

rowId2 <- 1

for (rowId2 in 1:nrow(dataSummary)) {
  
  plotData <- dplyr::filter(mergedData2, comparison == dataSummary$degId[rowId2]) %>% 
    dplyr::mutate(
      SM_CLUSTER = dplyr::if_else(
        condition = SM_ID %in% geneset$SM_ID[which(geneset$geneId == dataSummary$geneId[rowId2])],
        true = paste("*", SM_CLUSTER),
        false = SM_CLUSTER
      )
    ) %>% 
    dplyr::group_by(SM_CLUSTER) %>% 
    dplyr::arrange(chr, start, .by_group = TRUE) %>% 
    dplyr::mutate(index = row_number()) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(chr, start) %>% 
    dplyr::mutate(
      SM_CLUSTER = forcats::as_factor(SM_CLUSTER)
    )
  
  pltTitle <- paste(
    dataSummary$OESMTF_name[rowId2], " binding and log2(OE/WT) for all SM clusters", sep = ""
  )
  
  pt_bindingLfc <- ggplot(
    data = plotData,
  ) +
    geom_tile(
      mapping = aes(x = index, y = "DEG", fill = log2FoldChange, alpha = significance),
      size = 0.2, height = 1, color = "grey") +
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
      values = c("bound" = "black", "not-bound" = NA),
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
      facets = . ~ SM_CLUSTER, scales = "free_y",
      ncol = 6, dir = "v"
    ) +
    ggtitle(pltTitle) +
    guides(color = guide_legend(override.aes = list(shape = 22, size = 2, fill = "grey"))) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.spacing = unit(0.15, "lines"),
      strip.background = element_rect(fill="white", size = 0.2),
      strip.text.x = element_text(size = 13, hjust = 0),
      legend.text = element_text(size = 14),
      legend.position = "bottom",
      legend.title = element_text(size = 14, face = "bold"),
      plot.margin = unit(rep(0.2, 4), "cm")
    )
  
  print(pt_bindingLfc)
  
}

# ggsave(filename = paste(outPrefix, ".tiles_eg.png", sep = ""), plot = pt_bindingLfc,
#        width = 15, height = 8)

dev.off()


g <- ggplot_gtable(ggplot_build(pt_bindingLfc))

stript <- which(grepl('strip-t', g$layout$name))
fills <- c("red","green","blue","yellow")

