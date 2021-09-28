suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(RColorBrewer))


## generate TF summary plot heatmap
## TF ChIPseq correlation heatmap +
## TF peak count bar chart +
## peak occupancy heat map

rm(list = ls())

##################################################################################

analysisName <- "TF_binding"
outDir <- here::here("analysis", "09_TF_binding", "01_summary")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
file_peakSummary <- here::here("analysis", "02_QC_TF", "TF_ChIP_summary.best_replicates.tab")

TF_dataPath <- here::here("data", "TF_data")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

cutoff_macs2Pval <- 20

##################################################################################
productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::filter(has_polII_ChIP == "has_data", has_TF_ChIP == "has_data", copyNumber == "sCopy")

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = productionData$tfId,
  dataPath = TF_dataPath)

tfInfoList <- purrr::transpose(tfInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

tfSummary <- suppressMessages(readr::read_tsv(file = file_peakSummary))
# unname(AnnotationDbi::mapIds(x = orgDb, keys = anDf$SM_TF, column = "GENE_NAME", keytype = "GID"))

##################################################################################
## generate combinatorial binding matrix using
mat <- chipmine::combinatorial_binding_matrix(
  sampleInfo = tfInfo, peakFormat = "narrowPeak", summitRegion = 50
)

mat$width = mat$end - mat$start


## select only the regions where at least one sample has peak with -log10pvalue >= 20
matFiltered <- dplyr::filter_at(
  .tbl = mat, .vars = vars(starts_with("peakPval.")),
  .vars_predicate = any_vars(. >= cutoff_macs2Pval)
)

# glimpse(matFiltered)

ggplot(data = matFiltered, mapping = aes(x = width)) +
  geom_histogram(bins = 50) +
  labs(
    title = paste("Region width density plot, n =", nrow(matFiltered))
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, face = "bold")
  )

## peak occupancy binary matrix
matPeakOccupancy <- dplyr::select(matFiltered, name, starts_with("overlap.")) %>% 
  dplyr::mutate_if(is.logical, as.numeric) %>% 
  tibble::column_to_rownames(var = "name") %>% 
  as.matrix()

colnames(matPeakOccupancy) <- stringr::str_replace(
  string = colnames(matPeakOccupancy), pattern = "overlap.", replacement = ""
)


## transpose the matrices to plot data horizontally
matPeakOccupancy <- t(matPeakOccupancy)


## get the normalized coverage for peak regions
mergedRegions <- GenomicRanges::makeGRangesFromDataFrame(df = matFiltered)
mergedRegions$name <- matFiltered$name

regionCov <- chipmine::region_coverage_matrix(regions = mergedRegions, exptInfo = tfInfo)

## z-score of coverage across samples for each region
regionCovMat <- tibble::column_to_rownames(regionCov, var = "name") %>% 
  as.matrix() %>% 
  chipmine::scale_matrix_rows(add_attr = F)


## correlation matrix between all TF ChIPseq coverage
corrMat <- cor(x = regionCovMat)

anDf <- tibble::tibble(sampleId = rownames(corrMat)) %>% 
  dplyr::left_join(y = tfSummary, by = "sampleId")

anDf$geneName <- AnnotationDbi::mapIds(
  x = orgDb, keys = anDf$SM_TF, column = "GENE_NAME", keytype = "GID"
)

## row clustering for samples
rowClust <- hclust(dist(matPeakOccupancy))


## peak count bar chart annotation
an_peakCount <- HeatmapAnnotation(
  peakCount = anno_barplot(
    x = anDf$peaks_pval20,
    axis_param = list(gp = gpar(fontsize = 14))
  ),
  which = "row",
  width = unit(3, "cm"),
  annotation_name_side = "top",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontface = "bold", fontsize = 16),
  annotation_label = c("Peak count")
)


## TF ChIPseq correlation heatmap
ht_corr <- ComplexHeatmap::Heatmap(
  matrix = corrMat,
  col = colorRamp2(breaks = c(-1, 0, 1), colors = c("#276419", "#f7f7f7", "#8e0152")),
  cluster_rows = rowClust, cluster_columns = rowClust,
  show_column_dend = FALSE,
  column_title = "TF ChIPseq correlation",
  column_title_gp = gpar(fontface = "bold", fontsize = 16),
  width = 6,
  # row_names_gp = gpar(fontsize = 12),
  column_labels = anDf$geneName, column_names_gp = gpar(fontsize = 13),
  show_row_names = FALSE,
  # row_dend_reorder = T,
  # column_dend_reorder = T,
  heatmap_legend_param = list(
    title = "correlation", title_gp = gpar(fontsize = 14, fontface = "bold"),
    direction = "horizontal", title_position = "leftcenter",
    legend_width = unit(5, "cm"), grid_height = unit(1, "cm"),
    labels_gp = gpar(fontsize = 14)
  )
)

## peak occupancy heat map
ht_occupancy <- Heatmap(
  matrix = matPeakOccupancy,
  name = "matPeakOccupancy",
  col = c("1" = "#ff7f00", "0" = "black"),
  width = 7,
  column_title = "Peak occupancy",
  column_title_gp = gpar(fontface = "bold", fontsize = 16),
  # column_km = 6, column_gap = unit(0, "mm"),
  cluster_columns = TRUE, show_column_dend = FALSE,
  cluster_rows = FALSE,
  show_column_names = FALSE,
  row_labels = anDf$geneName, row_names_side = "left",
  row_names_gp = gpar(fontsize = 13),
  heatmap_legend_param = list(
    title = "",
    at = c(1), labels = "has peak",
    grid_height = unit(1, "cm"), grid_width = unit(1, "cm"),
    direction = "horizontal", title_position = "leftcenter",
    labels_gp = gpar(fontsize = 14)
  ),
  use_raster = TRUE, raster_quality = 5,
)


ptList <- ht_corr + ht_occupancy + an_peakCount

png(filename = paste(outPrefix, ".overall_comparison.png", sep = ""), width = 6500, height = 3500, res = 350)

draw(
  ptList,
  # padding = unit(c(2, 0.5, 0.5, 0.5), "cm"),
  heatmap_legend_side = "bottom",
  legend_gap = unit(10, "cm"),
  gap = unit(2, "mm"),
  auto_adjust = FALSE
)

dev.off()



