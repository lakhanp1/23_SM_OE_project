suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))

## polII DEG heatmap and significant DEG bar chart for SM clusters

rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################

analysisName <- "SM_cluster_DEG"
outDir <- here::here("analysis", "08_polII_diff_downstream", "SM_cluster_DEG")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

if(!dir.exists(outDir)){
  dir.create(outDir, recursive = T)
}

diffDataPath <- here::here("analysis", "07_polII_diff")
file_sampleInfo <- here::here("data", "reference_data", "polII_sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_degIds <- here::here("data", "reference_data", "polII_DEG_ids.txt")

file_clusterInfo <- here::here("data", "reference_data", "raw_data_summary.SM_cluster.tab")
file_backboneInfo <- here::here("data", "reference_data", "raw_data_summary.TF.tab")

useAllGroupsSamples <- FALSE

cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF
col_geneId <- "GID"

###########################################################################
dataSummary_SM <- suppressMessages(readr::read_tsv(file = file_clusterInfo)) %>% 
  dplyr::mutate(
    has_tf = forcats::fct_relevel(has_tf, "has_tf", "no_tf"),
    tf_data = forcats::fct_relevel(tf_data, "has_data", "no_data"),
    polII_data = forcats::fct_relevel(polII_data, "has_data", "no_data"),
  ) %>% 
  dplyr::arrange(has_tf, polII_data) %>% 
  dplyr::mutate(
    tf_with_chip = forcats::as_factor(tf_with_chip),
    tf_with_polII = forcats::as_factor(tf_with_polII)
  )

dataSummary_tf <- suppressMessages(readr::read_tsv(file = file_backboneInfo))

degIds <- suppressMessages(readr::read_tsv(file = file_degIds))

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% degIds$comparison)

rnaseqInfo <- dplyr::left_join(x = rnaseqInfo, y = dataSummary_tf, by = c("SM_TF" = "geneId")) %>% 
  dplyr::mutate(
    SM_cluster_tf = dplyr::if_else(
      condition = is.na(SM_ID), true = "no", false = "yes", missing = "no"
    ),
    SM_cluster_tf = forcats::fct_relevel(SM_cluster_tf, "yes", "no")
  ) %>% 
  dplyr::arrange(SM_cluster_tf)


geneInfo <- GenomicFeatures::genes(txDb) %>% 
  as.data.frame() %>% 
  dplyr::mutate(chr = as.character(seqnames)) %>% 
  dplyr::select(geneId = gene_id, chr, start, end, strand)

smGenes <- AnnotationDbi::select(
  x = orgDb, keys = keys(x = orgDb, keytype = "SM_CLUSTER"),
  columns = c("SM_ID", "GID"), keytype = "SM_CLUSTER"
) %>% 
  dplyr::rename(geneId = GID) %>% 
  tidyr::unite(SM_ID, geneId, col = "key", sep = ".", remove = FALSE) %>% 
  dplyr::left_join(y = dataSummary_SM, by = "SM_ID") %>% 
  dplyr::left_join(y = geneInfo, by = "geneId") %>% 
  dplyr::arrange(tf_with_polII, SM_ID, chr, start) %>% 
  dplyr::mutate(
    SM_CLUSTER = forcats::as_factor(SM_CLUSTER),
    SM_ID = forcats::as_factor(SM_ID)
  )


degData <- NULL

## get polII DEGs
row <- 1

for (row in 1:nrow(rnaseqInfo)) {
  tmpDf <- suppressMessages(readr::read_tsv(file = rnaseqInfo$deg[row])) %>% 
    dplyr::select(geneId, log2FoldChange, pvalue, padj) %>% 
    dplyr::mutate(comparison = rnaseqInfo$comparison[row])
  
  tmpDf <- dplyr::left_join(x = smGenes, y = tmpDf, by = "geneId")
  
  degData <- dplyr::bind_rows(degData, tmpDf)
}


degData <- degData %>% 
  dplyr::mutate(
    deg = dplyr::case_when(
      padj <= cutoff_fdr & log2FoldChange > 0 ~ "up",
      padj <= cutoff_fdr & log2FoldChange < 0 ~ "down",
      TRUE ~ "noDEG"
    ),
    deg = forcats::fct_relevel(deg, "up", "down", "noDEG"),
    comparison = forcats::fct_relevel(comparison, rnaseqInfo$comparison),
    log2FoldChange = if_else(condition = padj > cutoff_fdr, true = 0, false = log2FoldChange, missing = 0)
  )


lfcDf <- tidyr::pivot_wider(
  data = degData,
  id_cols = comparison,
  names_from = key,
  values_from = log2FoldChange,
  values_fill = list(key = 0)
)


lfcMat <- as.data.frame(lfcDf) %>% 
  tibble::column_to_rownames(var = "comparison") %>% 
  as.matrix()

if(!identical(x = colnames(lfcMat), y = smGenes$key)){
  stop("column annotation data.frame order does not match with heatmap matrix column order")
}

if(!identical(x = rownames(lfcMat), y = rnaseqInfo$comparison)){
  stop("row annotation data.frame order does not match with heatmap matrix row order")
}


an_top <- HeatmapAnnotation(
  polII_info = smGenes$tf_with_polII,
  which = "column",
  col = list(polII_info = c("has_tf:has_data" = "green", "has_tf:no_data" = "grey", "no_tf:no_data" = "black")),
  annotation_name_side = "left",
  annotation_legend_param = list(
    title = "has polII ChIPseq data",
    nrow = 1, title_position = "topcenter",
    grid_height = unit(1, "cm"), grid_width = unit(1, "cm"),
    labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14, fontface = "bold")
  )
)

ht <- ComplexHeatmap::Heatmap(
  matrix = lfcMat,
  name = "log2FoldChange",
  col = colorRamp2(breaks = c(-2, 0, 2), colors = c("#313695", "white", "#d73027")),
  row_title = "polII DEGs",
  column_title_side = "bottom",
  column_title_rot = 90,
  top_annotation = an_top,
  show_row_names = TRUE, row_names_side = "left",
  row_labels = stringr::str_replace(
    string = rownames(lfcMat), pattern = "_sCopy_OE_vs_MH11036", replacement = "_OE"
  ),
  row_dend_reorder = TRUE, column_dend_reorder = TRUE,
  show_column_names = FALSE, show_column_dend = FALSE,
  column_split = smGenes$SM_CLUSTER,
  row_split = rnaseqInfo$SM_cluster_tf,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.segments(x0 = x - width*0.5, y0 = y - height*0.5,
                  x1 = x + width*0.5, y1 = y - height*0.5,
                  gp = gpar(col = "grey")
    )
  },
  border = TRUE,
  cluster_column_slices = FALSE,
  heatmap_legend_param = list(
    title = "log2(fold-change)",
    direction = "horizontal", title_position = "topcenter",
    legend_width = unit(4, "cm"), grid_height = unit(1, "cm"), 
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 14)
  ),
  use_raster = FALSE
)


an_row <- HeatmapAnnotation(
  is_cluster_tf = rnaseqInfo$SM_cluster_tf,
  which = "row",
  col = list(is_cluster_tf = c("yes" = "green", "no" = "black")),
  annotation_legend_param = list(
    title = "SM cluster has TF",
    nrow = 1, title_position = "topcenter",
    grid_height = unit(1, "cm"), grid_width = unit(1, "cm"),
    labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14, fontface = "bold")
  )
)


plList <- ht + an_row

png(filename = paste(outPrefix, ".log2FC.heatmap.png", sep = ""), width = 8000, height = 5000, res = 400)

draw(
  plList,
  column_title = "polII ChIPseq DEGs log2FoldChange for SM cluster genes",
  column_title_gp = gpar(fontsize = 18, fontface = "bold"),
  heatmap_legend_side = "bottom",
  merge_legend = TRUE,
  padding = unit(c(10, 2, 2, 10), "mm")
)

dev.off()


###########################################################################


pt_bar <- degData %>%
  # dplyr::filter(SM_ID %in% c("cluster_01", "cluster_02", "cluster_03")) %>%
  # dplyr::filter(comparison %in% c("AN0148_sCopy_OE_vs_MH11036", "AN0153_sCopy_OE_vs_MH11036", "AN0533_sCopy_OE_vs_MH11036")) %>%
  ggplot(mapping = aes(x = comparison)) +
  geom_bar(mapping = aes(fill = deg), stat = "count",
           position = position_stack(reverse = TRUE)) +
  scale_fill_manual(
    values = c("up" = "#d73027", "noDEG" = "white", "down" = "#313695"),
    breaks = c("up", "down", "noDEG")
  ) +
  labs(
    title = "Significant DEG count bar plot in each SM cluster in polII OE/WT comparison"
  ) +
  scale_y_continuous(expand = expansion(add = 0)) +
  scale_x_discrete(
    labels = structure(
      stringr::str_replace(
        string = rnaseqInfo$comparison,
        pattern = "_sCopy_OE_vs_MH11036",
        replacement = "_polII"),
      names = rnaseqInfo$comparison
    ),
    expand = expansion(add = 0)
  ) +
  facet_grid(rows = vars(SM_CLUSTER), scales = "free") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    # panel.grid = element_blank(),
    panel.spacing = unit(0.5, "mm"),
    strip.background = element_rect(fill = "white"),
    strip.text.y = element_text(angle = 0, hjust = 0),
    strip.text.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"),
    legend.key = element_rect(color="black"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold")
  )


png(filename = paste(outPrefix, ".deg.barplot.png", sep = ""), width = 5000, height = 7000, res = 450)
# pdf(file = paste(outPrefix, ".deg.barplot.pdf", sep = ""), width = 12, height = 15)
print(pt_bar)
dev.off()



pt_barUp <- pt_bar +
  scale_fill_manual(
    values = c("up" = "#d73027", "noDEG" = "white", "down" = "white"),
    breaks = c("up", "down", "noDEG")
  )


png(filename = paste(outPrefix, ".deg_up.barplot.png", sep = ""), width = 5000, height = 7000, res = 450)
# pdf(file = paste(outPrefix, ".deg_up.barplot.pdf", sep = ""), width = 12, height = 15)
print(pt_barUp)
dev.off()








