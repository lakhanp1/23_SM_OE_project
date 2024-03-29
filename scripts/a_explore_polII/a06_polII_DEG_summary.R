suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))

## plot combined log2FoldChange heatmap of all SM OE polII strains

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################

analysisName <- "polII_DEGs_summary"
outDir <- here::here("analysis", "08_polII_analysis", "01_polII_DEGs_summary")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

if(!dir.exists(outDir)){
  dir.create(outDir, recursive = T)
}

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
diffDataPath <- here::here("analysis", "06_polII_diff")
file_sampleInfo <- here::here("data", "reference_data", "polII_sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_polIIDegs <- here::here("analysis", "06_polII_diff", "polII_DEGs.fpkm_filtered.combined.tab")

useAllGroupsSamples <- FALSE

cutoff_fpkm <- 10
cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1

orgDb <- org.Anidulans.FGSCA4.eg.db
col_geneId <- "GID"

col_lfc <- "log2FoldChange"
col_pval <- "padj"

###########################################################################
productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::filter(has_polII_ChIP == "has_data", has_TF_ChIP == "has_data", copyNumber == "sCopy")

smInfo <- AnnotationDbi::select(
  x = orgDb, keys = unique(na.exclude(productionData$SM_ID)),
  columns = c("SM_CLUSTER"), keytype = "SM_ID"
)

productionData <- dplyr::left_join(x = productionData, y = smInfo, by = "SM_ID")

degLabels <- structure(productionData$SMTF_name, names = productionData$degId)

degData <- suppressMessages(readr::read_tsv(file = file_polIIDegs)) %>% 
  dplyr::filter(comparison %in% productionData$degId)


degData <- dplyr::mutate(
  degData,
  diff = dplyr::case_when(
    !!sym(col_pval) <= cutoff_fdr & !!sym(col_lfc) <= cutoff_down ~ "down",
    !!sym(col_pval) <= cutoff_fdr & !!sym(col_lfc) >= cutoff_up ~ "up",
    TRUE ~ "noDEG"
  )
) %>% 
  dplyr::filter(diff != "noDEG", maxFpkm >= cutoff_fpkm)


## a dataframe with DEG status column for each comparison
diffDf <- tidyr::pivot_wider(
  data = degData,
  id_cols = c(geneId),
  names_from = comparison,
  values_from = diff,
  values_fill = list(diff = "noDEG")
)

readr::write_tsv(x = diffDf, file = paste(outPrefix, ".combined_DEGs.tab", sep = ""))

###########################################################################
## DEG stats bar plot
degStats <- dplyr::group_by(degData, comparison, diff) %>% 
  dplyr::count() %>% 
  tidyr::pivot_wider(
    names_from = diff, values_from = n, values_fill = list(n = 0)
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(total = up + down) %>% 
  dplyr::left_join(
    y = dplyr::select(productionData, degId, SMTF, SMTF_name, SM_CLUSTER),
    by = c("comparison" = "degId")
  )


readr::write_tsv(x = degStats, file = paste(outPrefix, ".stats.tab", sep = ""))


ptTheme <- theme_bw() +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 16),
    axis.title = element_blank(),
    title = element_text(size = 18),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = c(0.1, 0.9), legend.justification = c(0, 1),
    legend.text = element_text(size = 14)
  )


statsPlotData <- dplyr::mutate(degStats, down = -1*down) %>% 
  tidyr::pivot_longer(
    cols = c("up", "down"), names_to = "degGroup", values_to = "n"
  ) %>% 
  dplyr::arrange(desc(total)) %>% 
  dplyr::mutate(
    # geneLabel = forcats::as_factor(geneLabel),
    SMTF_name = forcats::as_factor(SMTF_name),
    labelHjust = if_else(condition = sign(n) == -1, true = 1, false = 0)
  )

pt_degStats <- ggplot(data = statsPlotData, mapping = aes(y = SMTF_name)) +
  geom_bar(mapping = aes(x = n, fill = degGroup), stat = "identity") +
  geom_text(
    mapping = aes(x = n+(20*sign(n)), label = abs(n), hjust = labelHjust)
  ) +
  geom_text(
    mapping = aes(x = max(n)+300, label = SM_CLUSTER), hjust = 0,
  ) +
  scale_fill_manual(
    values = c("up" = "#d73027", "down" = "#313695")
  ) +
  scale_x_continuous(
    breaks = seq(-1500, 1000, 500),
    expand = expansion(add = c(300, 1100))
  ) +
  labs(title = "TFOE/WT polII ChIPseq DEG count") +
  ptTheme

ggsave(filename = paste(outPrefix, ".stats.png", sep = ""), plot = pt_degStats, width = 12, height = 8)

###########################################################################
## heatmap of log2FoldChanges
lfcDf <- tidyr::pivot_wider(
  data = degData,
  id_cols = c(geneId),
  names_from = comparison,
  values_from = !!sym(col_lfc),
  values_fill = list(0) %>% purrr::set_names(nm = col_lfc)
)

readr::write_tsv(x = lfcDf, file = paste(outPrefix, ".log2FoldChange.tab", sep = ""))

lfcMat <- as.data.frame(lfcDf) %>% 
  tibble::column_to_rownames(var = "geneId") %>% 
  as.matrix()

colAnDf <- tibble::tibble(comparison = colnames(lfcMat)) %>% 
  dplyr::left_join(y = degStats, by = "comparison") %>% 
  dplyr::mutate(
    total = up + down,
    copy = stringr::str_replace(
      string = comparison, pattern = ".*_(\\w+Copy)_OE.*", replacement = "\\1"
    ),
    fraction = label_percent(accuracy = 0.01)(total/nrow(lfcMat))
  ) %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames(var = "comparison")


htAn <- ComplexHeatmap::HeatmapAnnotation(
  copy = colAnDf$copy,
  up = anno_text(
    x = colAnDf$up, rot = 90, just = "center", location = 0.5,
    gp = gpar(fill = "#ffb3c4", col = "black", border = "black"),
    height = max_text_width(text = colAnDf$up)*1.5
  ),
  down = anno_text(
    x = colAnDf$down, rot = 90, just = "center", location = 0.5,
    gp = gpar(fill = "#d8daf3", col = "black", border = "black"),
    height = max_text_width(text = colAnDf$down)*1.5
  ),
  fraction = anno_text(
    x = colAnDf$fraction, rot = 90, just = "center", location = 0.5,
    gp = gpar(col = "black", border = "black"),
    height = max_text_width(text = colAnDf$fraction)*1.2
  ),
  which = "column",
  col = list(
    copy = c("sCopy" = "#ff7f00", "msCopy" = "#8dd3c7")
  ),
  annotation_label = c(
    "copy number", "DEG Up", "DEG down", "fraction"
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "right",
  annotation_legend_param = list(
    copy = list(
      title = "Copy number", ncol = 1, title_position = "leftcenter-rot",
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    )
  )
)


ht_lfc <- ComplexHeatmap::Heatmap(
  matrix = lfcMat,
  name = "log2FoldChange",
  col = colorRamp2(breaks = c(-3, 0, 3), colors = c("#313695", "white", "#a50026")),
  top_annotation = htAn,
  row_title = "genes",
  column_title = "log2FoldChange heatmap of all significant DEGs in all OE/WT comparisons", 
  column_title_gp = gpar(fontsize = 18, fontface = "bold"),
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 16),
  row_dend_reorder = TRUE, column_dend_reorder = TRUE,
  column_labels = degLabels[colnames(lfcMat)],
  heatmap_legend_param = list(
    title = "log2(fold-change)",  title_position = "leftcenter-rot",
    legend_height = unit(4, "cm"), title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 14)
  ),
  use_raster = FALSE
)

htList <- ht_lfc

png(filename = paste(outPrefix, ".combined_lfc_heatmap.png", sep = ""), width = 6000, height = 4000, res = 400)
draw(
  htList,
  merge_legends = TRUE,
  legend_gap = unit(3, "cm"),
  padding = unit(c(0.5, 0.5, 0.5, 1), "cm")
)
dev.off()

###########################################################################


