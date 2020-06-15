suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(argparse))



## DESeq2 and RankComp2 DEG comparison

rm(list = ls())

##################################################################################

rankcompPath <- here::here("analysis", "07_polII_rank_diff")
diffDataPath <- here::here("analysis", "07_polII_diff")

file_sampleInfo <- here::here("data", "reference_data", "polII_sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_fpkm <- here::here("data", "polII_data", "polII_signal_matrix.tab")

cutoff_rank <- 5
cutoff_rankUp <- cutoff_rank
cutoff_rankDown <- cutoff_rank * -1
cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1

##################################################################################
## optionally run for multiple samples using command line arguments

parser <- ArgumentParser(
  description = "This script automates the DEG comparison by DESeq2 and RankComp2 algorithms."
)

parser$add_argument(
  "-d", "--deg", action="store",
  dest = "deg", type = "character", nargs = 1, required = TRUE,
  help = "DEG comparison ID. This value should be present in the column comparison of config file"
)

# parser$print_help()

# analysisName <- "AN0153_sCopy_OE_vs_MH11036"

args <- parser$parse_args()

analysisName <- args$deg

##################################################################################

outDir <- paste(rankcompPath, "/", analysisName, sep = "")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

diffInfo <- suppressMessages(readr::read_tsv(file = file_RNAseq_info)) %>% 
  dplyr::filter(comparison == analysisName)

samples <- unlist(stringr::str_split(string = diffInfo$samples, pattern = ";"))

sampleInfo <- read.table(file = file_sampleInfo, header = T, sep = "\t", stringsAsFactors = F)
fpkmMat <- suppressMessages(readr::read_tsv(file = file_fpkm))

##################################################################################


file_deseq <- paste(diffDataPath, "/", analysisName, "/", analysisName, ".DEG_all.txt", sep = "")
file_normCounts <- paste(diffDataPath, "/", analysisName, "/", analysisName, ".normCounts.tab", sep = "")
file_rankComp <- paste(rankcompPath, "/", analysisName, "/", analysisName, ".RankComp2.tab", sep = "")

normCounts <- suppressMessages(readr::read_tsv(file = file_normCounts))

deseqData <- suppressMessages(readr::read_tsv(file = file_deseq)) %>% 
  dplyr::select(geneId, log2FoldChange, deseq2_pval = pvalue, deseq2_padj = padj) %>% 
  dplyr::mutate(
    deseq2_diff = dplyr::case_when(
      deseq2_padj <= cutoff_fdr & log2FoldChange <= cutoff_down ~ "down",
      deseq2_padj <= cutoff_fdr & log2FoldChange >= cutoff_up ~ "up",
      TRUE ~ "noDEG"
    )
  ) %>% 
  dplyr::left_join(y = normCounts, by = "geneId")

rankcompData <- suppressMessages(readr::read_tsv(file = file_rankComp)) %>% 
  dplyr::rename(rankcomp_fdr = fdr, rankcomp_diff = diff)



mtSampleInfo <- dplyr::filter(sampleInfo, !!sym(diffInfo$design) %in% diffInfo$group1) %>% 
  dplyr::filter(sampleId %in% samples) %>% 
  dplyr::mutate(index = paste("MT", row_number(), sep = ""))

wtSampleInfo <- dplyr::filter(sampleInfo, !!sym(diffInfo$design) %in% diffInfo$group2) %>% 
  dplyr::filter(sampleId %in% samples) %>% 
  dplyr::mutate(index = paste("WT", row_number(), sep = ""))

exptInfo <- dplyr::bind_rows(wtSampleInfo, mtSampleInfo) %>% 
  dplyr::mutate(rankCol = paste("rank.", sampleId, sep = ""))

fpkms <- dplyr::select(fpkmMat, geneId, !!!samples) %>% 
  dplyr::rename(!!!structure(syms(samples), names = paste("fpkm.", samples, sep = "")))

rankcompData <- dplyr::left_join(x = rankcompData, y = fpkms, by = "geneId")

mergedData <- dplyr::left_join(x = deseqData, y = rankcompData, by = "geneId") %>% 
  dplyr::mutate(
    deseq2_diff = forcats::fct_relevel(deseq2_diff, "up", "noDEG", "down"),
    rankcomp_diff = forcats::fct_relevel(rankcomp_diff, "up", "noDEG", "down"),
    consensus = forcats::fct_relevel(consensus, "+", "+-", "-"),
    group = paste("DESeq2:", deseq2_diff, "; RankComp2:", rankcomp_diff,
                  "; rank change direction:", consensus, sep = "")
  ) %>% 
  dplyr::filter(deseq2_diff != "noDEG" | rankcomp_diff != "noDEG")

# glimpse(mergedData)
dplyr::group_by(mergedData, deseq2_diff, rankcomp_diff, consensus, group) %>% 
  dplyr::summarise(n = n())

readr::write_tsv(x = mergedData, path = paste(outPrefix, ".deseq2.rankcomp.comp.tab", sep = ""))

#########################
## prepare data

## deseq2 fold change
lfcMat <- dplyr::select(mergedData, geneId, log2FoldChange) %>% 
  tibble::column_to_rownames(var = "geneId") %>% 
  as.matrix()

## DEG status
diffMat <- dplyr::select(mergedData, geneId, deseq2_diff, rankcomp_diff, consensus) %>% 
  tibble::column_to_rownames(var = "geneId")

## deseqs normalized count z-score
normcountMat <- dplyr::select(mergedData, geneId, !!!samples) %>% 
  tibble::column_to_rownames(var = "geneId") %>% 
  as.matrix() %>% 
  chipmine::scale_matrix_rows(add_attr = FALSE)

## FPKM z-score
fpkmMat <- dplyr::select(mergedData, geneId, starts_with("fpkm.")) %>% 
  tibble::column_to_rownames(var = "geneId") %>% 
  as.matrix() %>% 
  chipmine::scale_matrix_rows(add_attr = FALSE)

## FPKM ranks
rankMat <- dplyr::select(mergedData, geneId, !!!paste("rank.", samples, sep = "")) %>% 
  tibble::column_to_rownames(var = "geneId") %>% 
  as.matrix()

## rank diff
rankDiffMat <- dplyr::select(mergedData, geneId, contains("_vs_")) %>% 
  tibble::column_to_rownames(var = "geneId") %>% 
  as.matrix() 

## mean rank, max rank, min rank
rankSummaryMat <- dplyr::select(mergedData, geneId, minRankDiff, maxRankDiff, avgRankDiff) %>% 
  tibble::column_to_rownames(var = "geneId") %>% 
  as.matrix()

#########################
## generate heatmaps

ptTitle <- paste(
  "DESeq2 (DEGs = ", length(which(deseqData$deseq2_diff != "noDEG")),
  ") and RankComp2 (DEGs = ", length(which(rankcompData$rankcomp_diff != "noDEG")),
  ") DEG comparison for ", analysisName, " (n = ", nrow(mergedData), ")",
  sep = "")

## DESeq2 log2FoldChange 
ht_lfc <- Heatmap(
  matrix = lfcMat,
  col = colorRamp2(breaks = c(-2, 0, 2), colors = c("#313695", "white", "#d73027")),
  name = "lfcMat",
  # column_title = "DESeq2 log2FoldChange",
  show_row_names = FALSE,
  cluster_columns = FALSE,
  row_split = diffMat,
  row_title = NULL,
  cluster_row_slices = FALSE,
  row_title_side = "left",
  row_title_rot = 0,
  width = unit(1.5, "cm")
)

## Group text zoom annotation
panel_fun = function(index, nm) {
  pushViewport(viewport())
  grid.rect()
  grid.text(
    x = unit(0.98, "npc"), y = unit(0.5, "npc"),
    label = paste(
      "(n = ", length(index),") ",
      stringr::str_replace_all(
        string = unique(mergedData$group[index]), pattern = "; ", replacement = "\n"
      ),
      sep = ""
    ),
    hjust = 1, gp = gpar(fontsize = 8)
  )
  popViewport()
}

an_group <- anno_zoom(
  align_to = mergedData$group,
  panel_fun = panel_fun,
  which = "row", side = "left",
  size = 1,
  gap = unit(0.5, "cm"),
  width = unit(5, "cm"),
  link_width = unit(1.5, "cm"), link_height = unit(1.5, "cm"),
  link_gp = gpar(fill = "#bdbdbd")
)

an_right <- rowAnnotation(group = an_group)

## gene group annotation
an_diff <- HeatmapAnnotation(
  df = diffMat,
  which = "row",
  col = list(
    deseq2_diff = c("down" = "#313695", "noDEG" = "grey", "up" = "#d73027"),
    rankcomp_diff = c("down" = "#313695", "noDEG" = "grey", "up" = "#d73027"),
    consensus = c("-" = "#ff7f00", "+" = "#984ea3", "+-" = "#999999")
  ),
  simple_anno_size = unit(0.7, "cm")
)

## DESeq2 normalized read count z-score Heatmap
ht_normCount <- Heatmap(
  matrix = normcountMat,
  name = "normcountMat",
  col = colorRamp2(breaks = c(-2, 0, 2), colors = c("#8c510a", "#f5f5f5", "#01665e")),
  column_title = "z-score(DESeq2 count)",
  column_labels = stringr::str_replace(
    # string = colnames(normcountMat), pattern = "(.*)_16h_polII_ChIP(Mix.*)", replacement = "\\1_\\2"
    string = colnames(normcountMat), pattern = "16h_polII", replacement = "16h\npolII"
  ),
  show_row_names = FALSE,
  cluster_columns = FALSE,
  width = unit(5, "cm")
)

## FPKM z-score heatmap
ht_fpkm <- Heatmap(
  matrix = fpkmMat,
  name = "fpkmMat",
  col = colorRamp2(breaks = c(-2, 0, 2), colors = c("#8e0152", "#f7f7f7", "#276419")),
  column_title = "z-score(FPKM)",
  column_labels = stringr::str_replace(
    string = colnames(fpkmMat), pattern = "16h_polII", replacement = "16h\npolII"
  ),
  show_row_names = FALSE,
  cluster_columns = FALSE,
  width = unit(5, "cm")
)

## rank(FPKM) Heatmap
ht_rank <- Heatmap(
  matrix = rankMat,
  name = "rankMat",
  col = colorRamp2(
    # breaks = c(1000, 3000, 5000, 7000, 8000, 9000, 9500, 10000, 10500),
    breaks = c(1000, 2000, 3000, 4000, 5000, 6000, 7500, 9000, 10500),
    colors = RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")
  ),
  column_title = "rank(FPKM)",
  column_labels = stringr::str_replace(
    string = colnames(rankMat), pattern = "16h_polII", replacement = "16h\npolII"
  ),
  show_row_names = FALSE,
  cluster_columns = FALSE,
  width = unit(5, "cm")
)

## rank difference Heatmap
ht_rankDiff <- Heatmap(
  matrix = rankDiffMat,
  name = "rankDiffMat",
  col = colorRamp2(
    breaks = seq(-50, 50, length.out = 11),
    colors = RColorBrewer::brewer.pal(n = 11, name = "PiYG"),
  ),
  column_title = "rank difference",
  show_row_names = FALSE,
  cluster_columns = FALSE,
  width = unit(5, "cm")
)

## rank summary Heatmap
ht_rankSummary <- Heatmap(
  matrix = rankSummaryMat,
  name = "rankSummaryMat",
  col = colorRamp2(
    breaks = seq(-50, 50, length.out = 11),
    colors = RColorBrewer::brewer.pal(n = 11, name = "PiYG"),
  ),
  column_title = "rank diff summary",
  show_row_names = FALSE,
  cluster_columns = FALSE,
  width = unit(5, "cm")
)

htLits <- an_right + an_diff + ht_lfc + ht_normCount + ht_fpkm + ht_rank + ht_rankDiff + ht_rankSummary


png(filename = paste(outPrefix, ".deseq2.rankcomp.comp.png", sep = ""), width = 7000, height = 5000, res = 450)

draw(
  htLits,
  column_title = ptTitle,
  column_title_gp = gpar(fontsize = 16),
  show_row_dend = FALSE,
  row_dend_side = "left",
  row_sub_title_side = "left",
  row_gap = unit(2, "mm")
)

dev.off()


