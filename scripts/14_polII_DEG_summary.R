suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))

## plot combined log2FoldChange heatmap of all SM OE polII strains

rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################

analysisName <- "combined_polII_DEGs"
outDir <- here::here("analysis", "08_polII_diff_downstream", "combined_polII_DEGs")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

if(!dir.exists(outDir)){
  dir.create(outDir, recursive = T)
}

diffDataPath <- here::here("analysis", "06_polII_diff")
file_sampleInfo <- here::here("data", "reference_data", "polII_sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_degIds <- here::here("data", "reference_data", "polII_DEG_ids.txt")

useAllGroupsSamples <- FALSE

cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1

orgDb <- org.Anidulans.FGSCA4.eg.db
col_geneId <- "GID"

###########################################################################
degIds <- suppressMessages(readr::read_tsv(file = file_degIds))

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath)


rowId <- 1
degData <- NULL

for (rowId in 1:nrow(rnaseqInfo)) {
  
  tmpDf <- suppressMessages(readr::read_tsv(file = rnaseqInfo$deg[rowId])) %>% 
    dplyr::select(geneId, contrast, log2FoldChange, padj, GENE_NAME) %>% 
    dplyr::mutate(comparison = rnaseqInfo$comparison[rowId])
  
  degData <- dplyr::bind_rows(degData, tmpDf)
}


degData <- dplyr::mutate(
  degData,
  # contrast = stringr::str_replace(string = contrast, pattern = "(condition_|SM_TF_)", replacement = ""),
  diff = dplyr::case_when(
    padj <= cutoff_fdr & log2FoldChange <= cutoff_down ~ "down",
    padj <= cutoff_fdr & log2FoldChange >= cutoff_up ~ "up",
    TRUE ~ "noDEG"
  )
) %>% 
  dplyr::filter(diff != "noDEG")

degStats <- dplyr::group_by(degData, comparison, diff) %>% 
  dplyr::count() %>% 
  tidyr::pivot_wider(
    names_from = diff, values_from = n, values_fill = list(n = 0)
  )

readr::write_tsv(x = degStats, path = paste(diffDataPath, "/polII_DEG.stats.tab", sep = ""))

## a dataframe with DEG status column for each comparison
diffDf <- tidyr::pivot_wider(
  data = degData,
  id_cols = c(geneId),
  names_from = comparison,
  values_from = diff,
  values_fill = list(diff = "noDEG")
)

readr::write_tsv(x = diffDf, path = paste(outPrefix, ".diff.tab", sep = ""))


###########################################################################
## heatmap of log2FoldChanges
lfcDf <- tidyr::pivot_wider(
  data = degData,
  id_cols = c(geneId),
  names_from = comparison,
  values_from = log2FoldChange,
  values_fill = list(log2FoldChange = 0)
)

readr::write_tsv(x = lfcDf, path = paste(outPrefix, ".log2FoldChange.tab", sep = ""))

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
  row_dend_reorder = TRUE, column_dend_reorder = TRUE,
  column_labels = stringr::str_replace(
    string = colnames(lfcMat), pattern = "_vs_MH11036", replacement = ""
  ),
  heatmap_legend_param = list(
    title = "log2(fold-change)",  title_position = "leftcenter-rot",
    legend_height = unit(4, "cm"), title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 14)
  ),
  use_raster = FALSE
)

htList <- ht_lfc

png(filename = paste(outPrefix, ".lfc_heatmap.png", sep = ""), width = 6000, height = 4000, res = 400)
draw(
  htList,
  merge_legends = TRUE,
  legend_gap = unit(3, "cm"),
  padding = unit(c(0.5, 0.5, 0.5, 1), "cm")
)
dev.off()

###########################################################################

# keytypes(orgDb)

smGenes <- AnnotationDbi::select(
  x = orgDb, keys = keys(x = orgDb, keytype = "SM_GENE"),
  columns = c("SM_GENE", "SM_CLUSTER")
)

tfGenes <- AnnotationDbi::select(
  x = orgDb, keys = keys(x = orgDb, keytype = "TF_GENE"),
  columns = c("GID")
)

geneSets <- list(
  deg_up = unique(degData$geneId[degData$diff == "up"]),
  deg_down = unique(degData$geneId[degData$diff == "down"]),
  SM_genes = unique(smGenes$GID),
  tf = unique(tfGenes$GID)
)

# ## Venn diagram
# VennDiagram::venn.diagram(
#   x = geneSets,
#   filename = paste(outPrefix, ".SM_DEG_Venn.png"),
#   fill = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072"),
#   fontface = "bold", cex = 1.2, 
#   cat.fontface = "bold", cat.cex = 1.2
# )


## Upset plot
cm = make_comb_mat(geneSets)


set_size(cm)
comb_name(cm)
comb_size(cm)
comb_degree(cm)

pt <- UpSet(
  m = cm,
  pt_size = unit(7, "mm"), lwd = 3,
  comb_col = RColorBrewer::brewer.pal(n = 4, name = "Dark2")[comb_degree(cm)],
  set_order = names(geneSets),
  top_annotation = HeatmapAnnotation(
    foo = anno_empty(border = FALSE),
    "combSize" = anno_text(
      x = paste("(", comb_size(cm), ")", sep = ""),
      just = "left", location = 0.5, rot = 90
    ),
    "Intersection\nsize" = anno_barplot(
      x = comb_size(cm),
      border = FALSE,
      gp = gpar(fill = "black"),
      height = unit(2, "cm")
    ),
    annotation_name_side = "left",
    annotation_name_rot = 0
  ),
  right_annotation = upset_right_annotation(
    m = cm, bar_width = 0.5
  ),
  width = unit(15, "cm"), height = unit(6, "cm")
)


png(filename = paste(outPrefix, ".SM_DEG_Upset.png", sep = ""), width = 6000, height = 4000, res = 400)

draw(
  pt,
  column_title = "SM OE strain DEG overlap with SM cluster genes and TF",
)

dev.off()






