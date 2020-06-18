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

analysisName <- "TF_peak_summary"
outDir <- here::here("analysis", "09_TF_binding")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_tfSamples <- here::here("data", "reference_data", "tf_sample_ids.best_rep.txt")
file_peakSummary <- here::here("analysis", "02_QC_TF", "TF_ChIP_summary.best_replicates.tab")
TF_dataPath <- here::here("data", "TF_data")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

cutoff_macs2Pval <- 20

##################################################################################


tfSampleList <- suppressMessages(readr::read_tsv(file = file_tfSamples))

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSampleList$sampleId,
  dataPath = TF_dataPath)

# glimpse(tfInfo)

tfInfoList <- purrr::transpose(tfInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

tfSummary <- suppressMessages(readr::read_tsv(file = file_peakSummary))

##################################################################################
## generate combinatorial binding matrix using
mat <- combinatorial_binding_matrix(
  sampleInfo = tfInfo, peakFormat = "narrowPeak", summitRegion = 50
)

mat$width = mat$end - mat$start


## select only the regions where at least one sample has peak with -log10pvalue >= 20
matFiltered <- dplyr::filter_at(
  .tbl = mat, .vars = vars(starts_with("peakPval.")),
  .vars_predicate = any_vars(. >= cutoff_macs2Pval)
)

# glimpse(matFiltered)


## peak occupancy binary matrix
matPeakOccupancy <- dplyr::select(matFiltered, name, starts_with("overlap.")) %>% 
  dplyr::mutate_if(is.logical, as.numeric) %>% 
  tibble::column_to_rownames(var = "name") %>% 
  as.matrix()

colnames(matPeakOccupancy) <- stringr::str_replace(
  string = colnames(matPeakOccupancy), pattern = "overlap.", replacement = ""
)


## peakEnrichment matrix
enrichmentMat <- dplyr::select(matFiltered, name, starts_with("peakEnrichment")) %>% 
  dplyr::mutate_at(.vars = vars(starts_with("peakEnrichment.")), .funs = ~ replace_na(., 0)) %>% 
  tibble::column_to_rownames(var = "name") %>% 
  as.matrix() %>% 
  scale_matrix_columns(add_attr = F)

colnames(enrichmentMat) <- gsub(pattern = "peakEnrichment.", replacement = "", x = colnames(enrichmentMat))

## transpose the matrices to plot data horizontally
matPeakOccupancy <- t(matPeakOccupancy)


## get the normalized coverage for peak regions
mergedRegions <- GenomicRanges::makeGRangesFromDataFrame(df = matFiltered)
mergedRegions$name <- matFiltered$name

regionCov <- region_coverage_matrix(regions = mergedRegions, exptInfo = tfInfo)

## z-score of coverage across samples for each region
regionCovMat <- tibble::column_to_rownames(regionCov, var = "name") %>% 
  as.matrix() %>% 
  chipmine::scale_matrix_rows(add_attr = F)


## correlation matrix between all TF ChIPseq coverage
corrMat <- cor(x = regionCovMat)

anDf <- tibble::tibble(sampleId = rownames(corrMat)) %>% 
  dplyr::left_join(y = tfSummary, by = "sampleId")

## row clustering for samples
rowClust <- hclust(dist(matPeakOccupancy))


## peak count bar chart annotation
an_peakCount <- HeatmapAnnotation(
  peakCount = anno_barplot(
    x = anDf$peaks_pval20,
    axis_param = list(gp = gpar(fontsize = 14))
  ),
  which = "row",
  width = unit(5, "cm"),
  annotation_name_side = "top",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontface = "bold", fontsize = 16),
  annotation_label = c("Peak count\n")
)


## TF ChIPseq correlation heatmap
ht_corr <- ComplexHeatmap::Heatmap(
  matrix = corrMat,
  col = colorRamp2(breaks = c(-1, 0, 1), colors = c("#276419", "#f7f7f7", "#8e0152")),
  # right_annotation = an_peakCount,
  cluster_rows = rowClust, cluster_columns = rowClust,
  column_title = "TF ChIPseq correlation",
  column_title_gp = gpar(fontface = "bold", fontsize = 16),
  width = 6,
  # row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  row_labels = stringr::str_replace(
    string = rownames(corrMat), pattern = "_sCopy_OE_16h_(HA|HIS|FLAG)_.*", replacement = ""
  ),
  column_labels = stringr::str_replace(
    string = rownames(corrMat), pattern = "_sCopy_OE_16h_(HA|HIS|FLAG)_.*", replacement = ""
  ),
  row_dend_reorder = T,
  column_dend_reorder = T,
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
  col = c("1" = "#ff7f00", "0" = "#e6e6e6"),
  width = 5,
  column_title = "Peak occupancy",
  column_title_gp = gpar(fontface = "bold", fontsize = 16),
  column_km = 6, column_gap = unit(0, "mm"),
  cluster_columns = TRUE, show_column_dend = FALSE,
  cluster_rows = FALSE,
  show_column_names = FALSE,
  row_labels = stringr::str_replace(
    string = rownames(matPeakOccupancy), pattern = "_sCopy_OE_16h_(HA|HIS|FLAG)_.*", replacement = ""
  ),
  heatmap_legend_param = list(
    title = "Peak occupancy", title_gp = gpar(fontsize = 14),
    at = c(1), labels = "has peak",
    grid_height = unit(1, "cm"), grid_width = unit(1, "cm"),
    direction = "horizontal", title_position = "leftcenter"
  ),
  use_raster = TRUE, raster_quality = 5,
)


ptList <- ht_corr + an_peakCount + ht_occupancy

png(filename = paste(outPrefix, ".binding_matrix.png", sep = ""), width = 6000, height = 3500, res = 350)

draw(
  ptList,
  padding = unit(c(2, 0.5, 0.5, 0.5), "cm"),
  heatmap_legend_side = "bottom",
  legend_gap = unit(10, "cm"),
  gap = unit(2, "mm")
)

dev.off()






