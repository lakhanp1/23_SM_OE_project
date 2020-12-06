suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))

## this script generate figures showing polII OE/WT for all SM clusters for a TFOE data

rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################

analysisName <- "cross_cluster_DEG"
outDir <- here::here("analysis", "08_polII_diff_downstream", "03_cross_cluster_DEG")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

if(!dir.exists(outDir)){
  dir.create(outDir, recursive = T)
}

diffDataPath <- here::here("analysis", "06_polII_diff")
file_sampleInfo <- here::here("data", "reference_data", "polII_sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_degIds <- here::here("data", "reference_data", "production_data.polII_DEG_ids.txt")

useAllGroupsSamples <- FALSE

cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF
col_geneId <- "GID"

col_lfc <- "log2FoldChange"
col_pval <- "pvalue"

###########################################################################
degIds <- suppressMessages(readr::read_tsv(file = file_degIds))

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% degIds$degId)

clusterInfo <- AnnotationDbi::select(
  x = orgDb, keys = rnaseqInfo$SM_TF, columns = c("GENE_NAME", "SM_CLUSTER"), keytype = "GID"
) %>% 
  dplyr::group_by(GID) %>% 
  dplyr::summarise(
    geneName = unique(GENE_NAME),
    SM_cluster = paste(na.exclude(unique(SM_CLUSTER)), collapse = "/")
  ) %>% 
  dplyr::arrange(desc(SM_cluster))

rnaseqInfo <- dplyr::left_join(x = rnaseqInfo, y = clusterInfo, by = c("SM_TF" = "GID")) %>% 
  dplyr::mutate(
    geneLabel = paste(SM_TF, " (", geneName, ")", sep = ""),
    geneLabel = if_else(condition = SM_TF == geneName, true = SM_TF, false = geneLabel),
    degLabel = paste(geneLabel, ": ", SM_cluster, sep = "")
  ) %>% 
  dplyr::arrange(desc(SM_cluster))


## extract genes belonging to "GO:0022900 electron transport chain"
geneset <- AnnotationDbi::select(
  x = orgDb, keys = keys(x = orgDb, keytype = "SM_ID"),
  columns = c("GID", "GENE_NAME", "SM_CLUSTER"), keytype = "SM_ID"
) %>% 
  dplyr::rename(geneId = GID) %>% 
  dplyr::mutate(
    geneLabel = paste(geneId, "(", GENE_NAME, ")", sep = ""),
    geneLabel = if_else(condition = geneId == GENE_NAME, true = geneId, false = geneLabel)
  )

genePos <- as.data.frame(genes(txDb)) %>% 
  dplyr::select(geneId = gene_id, chr = seqnames, start, end, strand)

geneset <- dplyr::left_join(x = geneset, y = genePos, by = "geneId") %>% 
  dplyr::arrange(SM_ID, chr, start)


###########################################################################
## get DEG data for each SM_TF cluster from its own polII_DEG set
rowId <- 1
degData <- NULL


for (rowId in 1:nrow(rnaseqInfo)) {
  
  tmpDf <- suppressMessages(readr::read_tsv(file = rnaseqInfo$deg[rowId])) %>% 
    dplyr::select(geneId, !!col_lfc, shrinkLog2FC, pvalue, padj) %>% 
    dplyr::mutate(
      comparison = rnaseqInfo$comparison[rowId],
      SM_TF = rnaseqInfo$SM_TF[rowId],
      degLabel = rnaseqInfo$degLabel[rowId]
    )
  
  subData <- dplyr::left_join(
    x = geneset, y = tmpDf, by = "geneId"
  )
  
  degData <- dplyr::bind_rows(degData, subData)
  
}


# readr::write_tsv(x = degData, file = paste(outPrefix, ".combined_DEGs.tab", sep = ""))



###########################################################################
## combined plot

## for diagonal symmetry of polII DEG data and SM clusters
clusterLevels <- unname(
  AnnotationDbi::mapIds(
    x = orgDb, keys = rnaseqInfo$SM_TF, column = "SM_CLUSTER", keytype = "GID"
  )
)

clusterLevels <- unique(clusterLevels[!is.na(clusterLevels)])

smDegSummary <- dplyr::mutate(
  degData,
  diff = dplyr::case_when(
    !!sym(col_pval) <= cutoff_fdr & !!sym(col_lfc) <= cutoff_down ~ "down",
    !!sym(col_pval) <= cutoff_fdr & !!sym(col_lfc) >= cutoff_up ~ "up",
    TRUE ~ "noDEG"
  )
) %>% 
  dplyr::count(comparison, degLabel, SM_CLUSTER, diff) %>% 
  dplyr::add_count(comparison, degLabel, SM_CLUSTER, wt = n, name = "total") %>% 
  dplyr::mutate(
    fraction = n/total,
    diff = forcats::fct_relevel(diff, "up", "noDEG", "down"),
    degLabel = forcats::fct_relevel(degLabel, !!!rnaseqInfo$degLabel),
    SM_CLUSTER = forcats::fct_relevel(SM_CLUSTER, !!!clusterLevels)
  )


trim_label <- function(char){
  char <- stringr::str_sub(string = char, end = 15)
  char
}

# dplyr::filter(smDegSummary, comparison %in% c("AN0148_sCopy_OE_vs_WT", "AN0153_sCopy_OE_vs_WT")) %>%
#   dplyr::filter(SM_CLUSTER %in% c("No PKS/NRPS backbone 1", "AN0016 cluster",
#                                   "AN10297 cluster", "AN11065 cluster")) %>%

pt_degPi <- smDegSummary %>%
  ggplot() +
  coord_fixed() +
  geom_arc_bar(
    mapping = aes(x0 = 0, y0 = 0, r0 = 0, r = 1, amount = n, fill = diff),
    stat = "pie", linetype = "solid") +
  scale_fill_manual(
    name = "DEG group",
    values = c("up" = "#d73027", "down" = "#313695", "noDEG" = "grey")
  ) +
  labs(title = "SMTF OE DEG proportion for all SM clusters") +
  facet_grid(
    facets = degLabel ~ SM_CLUSTER,
    labeller = labeller(degLabel = trim_label, SM_CLUSTER = trim_label)
  ) +
  theme_no_axes() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.spacing = unit(0.05, "lines"),
    strip.background = element_rect(fill="white", size = 0.2),
    strip.text.x = element_text(size = 15, angle = 90, hjust = 0),
    strip.text.y = element_text(size = 15, angle = 0, hjust = 0),
    legend.text = element_text(size = 15),
    legend.position = "right",
    legend.title = element_text(size = 15, face = "bold", angle = 90)
  )


ggsave(filename = paste(outPrefix, ".pie.png", sep = ""), plot = pt_degPi, width = 20, height = 12)

###########################################################################

pdf(file = paste(outPrefix, ".TFOE_wise.pdf", sep = ""), width = 15, height = 8, onefile = TRUE)

rowId2 <- 1

for (rowId2 in 1:nrow(rnaseqInfo)) {
  
  plotData <- dplyr::filter(degData, comparison == rnaseqInfo$comparison[rowId2]) %>% 
    dplyr::group_by(SM_CLUSTER) %>% 
    dplyr::arrange(chr, start, .by_group = TRUE) %>% 
    dplyr::mutate(index = row_number()) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      significance = if_else(
        condition = !!sym(col_pval) <= cutoff_fdr, true = "significant",
        false = "non-significant", missing = "non-significant"
      )
    )
  
  pltTitle <- paste(rnaseqInfo$degLabel[rowId2], " || log2(", rnaseqInfo$SM_TF[rowId2],
                    "-OE/WT) for all SM clusters", sep = "")
  
  pt_tiles <- ggplot() +
    geom_tile(
      data = plotData,
      mapping = aes(x = index, y = SM_CLUSTER, fill = shrinkLog2FC, color = significance),
      size = 0.2, height = 1) +
    scale_fill_gradient2(
      name = paste("log2(", "OE/WT", ")", sep = ""),
      low = "#313695", mid = "#F7F7F7", high = "#d73027", midpoint = 0,
      limit = c(-2.5, 2.5), oob = scales::squish
    ) +
    scale_colour_manual(
      name = "p-value <= 0.05",
      values = c("significant" = "black", "non-significant" = "white")
    ) +
    scale_x_continuous(expand = expansion(add = c(0.0, 0.0))) +
    facet_wrap(
      facets = . ~ SM_CLUSTER, scales = "free_y",
      ncol = 5, dir = "v"
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
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),
      plot.margin = unit(rep(0.2, 4), "cm")
    )
  
  print(pt_tiles)
  
  # ggsave(filename = paste(outPrefix, ".tmp.png", sep = ""), plot = pt_tiles, width = 14, height = 8)
  
}

dev.off()

