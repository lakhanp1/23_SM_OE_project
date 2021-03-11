suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))


## this script:
## 1) use diffbind regions and generate profile plot around peak

rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################
outDir <- here::here("analysis", "11_aflR_AN7820_analysis", "a05_diffbind_DEG_integration")

analysisName <- "AflR_diffbind_vs_DEG"
outPrefix <- paste(outDir, "/", analysisName, sep = "")

## the denominator or WT in log2(fold_change) should be second
col_compare <- "diffbind_condition"
diffbindCompare <- c("low_Xylose", "high_Xylose")
degId <- "AN7820_sCopy_OE_vs_WT"


file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_diffbindRes <- here::here(
  "analysis", "11_aflR_AN7820_analysis", "a04_diffbind", "AflR_diffbind.annotation.filtered.tab"
)

TF_dataPath <- here::here("data", "TF_data")
diffDataPath <- here::here("analysis", "06_polII_diff")

orgDb <- org.Anidulans.FGSCA4.eg.db

cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1
##################################################################################
geneDesc <- suppressMessages(
  AnnotationDbi::select(x = orgDb, keys = keys(x = orgDb),
                        columns = c("DESCRIPTION"), keytype = "GID")
) %>% 
  dplyr::rename(geneId = GID)

smInfo <- suppressMessages(
  AnnotationDbi::select(x = orgDb, keys = keys(orgDb, keytype = "SM_CLUSTER"),
                        columns = c("GID", "SM_ID"), keytype = "SM_CLUSTER")) %>% 
  dplyr::group_by(GID) %>% 
  dplyr::mutate(SM_CLUSTER = paste(SM_CLUSTER, collapse = ";"),
                SM_ID = paste(SM_ID, collapse = ";")) %>% 
  dplyr::slice(1L) %>% 
  dplyr::ungroup()

geneInfo <- dplyr::left_join(x = geneDesc, y = smInfo, by = c("geneId" = "GID"))

glimpse(geneInfo)

##################################################################################
if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

diffbindRes <- suppressMessages(readr::read_tsv(file = file_diffbindRes)) %>% 
  dplyr::filter(pvalGood.all != 0) %>% 
  dplyr::filter(diffBind != "noDiff" & !is.na(geneId)) %>%
  dplyr::select(
    diffbind_region = name,
    diffbind_lfc = Fold,
    geneId,
    diffbind_FDR = FDR,
    diffbind_pval = p.value,
    diffbind_diff = diffBind
  )

deseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% degId)

normCounts <- suppressMessages(readr::read_tsv(file = deseqInfo$normCount))

deseqData <- suppressMessages(readr::read_tsv(file = deseqInfo$deseq2)) %>% 
  dplyr::select(geneId, deseq2_lfc = log2FoldChange, deseq2_pval = pvalue, deseq2_padj = padj) %>% 
  dplyr::mutate(
    deseq2_diff = dplyr::case_when(
      deseq2_padj <= cutoff_fdr & deseq2_lfc <= cutoff_down ~ "down",
      deseq2_padj <= cutoff_fdr & deseq2_lfc >= cutoff_up ~ "up",
      TRUE ~ "noDEG"
    )
  ) %>% 
  dplyr::left_join(y = normCounts, by = "geneId")


mergedData <- dplyr::left_join(x = diffbindRes, y = deseqData, by = "geneId") %>% 
  dplyr::mutate(
    deseq2_diff = forcats::fct_relevel(deseq2_diff, "up", "noDEG", "down"),
    diffbind_diff = forcats::fct_relevel(diffbind_diff, "up", "down")
  ) %>% 
  dplyr::left_join(y = geneInfo, by = "geneId") %>% 
  tidyr::unite(col = "rowId", diffbind_region, geneId, remove = F) %>% 
  dplyr::filter(!is.na(SM_ID))



diffbind_lfcMat <- dplyr::select(mergedData, rowId, diffbind_lfc) %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames(var = "rowId") %>% 
  as.matrix()

polII_lfcMat <- dplyr::select(mergedData, rowId, deseq2_lfc) %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames(var = "rowId") %>% 
  as.matrix()

polII_pval <- dplyr::mutate(
  mergedData,
  deseq2_significant = dplyr::if_else(
    condition = deseq2_padj < cutoff_fdr, true = "yes", false = "no", missing = "no"
  )
) %>% 
  dplyr::select(rowId, deseq2_significant)

##################################################################################
## DESeq2 log2FoldChange 
ht_diffbind <- Heatmap(
  matrix = diffbind_lfcMat,
  col = colorRamp2(breaks = c(-3, 0, 3), colors = c("#313695", "white", "#a50026")),
  name = "DiffBind",
  column_title = "DiffBind fold-change\n0.05% Xyl vs 1% Xylose",
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  row_split = mergedData$SM_CLUSTER,
  # row_title = NULL,
  # cluster_row_slices = FALSE,
  row_title_side = "left",
  row_title_rot = 0,
  width = unit(5, "cm")
)


ht_polII <- Heatmap(
  matrix = polII_lfcMat,
  col = colorRamp2(breaks = c(-3, 0, 3), colors = c("#313695", "white", "#a50026")),
  name = "polII ChIP",
  column_title = "polII ChIPseq fold-change\nAflR-OE vs WT",
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  row_split = mergedData$diffbind_diff,
  row_title = NULL,
  cluster_row_slices = FALSE,
  row_title_side = "left",
  row_title_rot = 0,
  width = unit(5, "cm")
)

htAn <- ComplexHeatmap::HeatmapAnnotation(
  significant = polII_pval$deseq2_significant,
  which = "row",
  col = list(
    significant = c("yes" = "black", "no" = "white")
  ),
  annotation_label = c("p-value"),
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  simple_anno_size = unit(2, "cm"),
  annotation_legend_param = list(
    significant = list(
      title = "p-value <= 0.05", ncol = 1,
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    )
  )
)

htList <- ht_diffbind + ht_polII + htAn

png(filename = paste(outPrefix, ".SM_lfc_heatmap.png", sep = ""), width = 5000, height = 3000, res = 400)
draw(
  htList,
  column_title = "AflR binding difference and polII ChIPseq comparison",
  column_title_gp = gpar(fontsize = 16),
  merge_legends = TRUE,
  padding = unit(c(0.5, 0.5, 0.5, 1), "cm")
)
dev.off()






