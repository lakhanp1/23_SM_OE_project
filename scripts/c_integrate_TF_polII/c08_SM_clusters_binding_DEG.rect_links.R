suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(gginnards))
suppressPackageStartupMessages(library(bezier))

## this script plots binding and polII-DEG data for all SMTFs over 71 SM clusters
## as rectangular link plots


rm(list = ls())

source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

##################################################################################

analysisName <- "SM_clusters_binding_DEG.curve_links"
outDir <- here::here("analysis", "10_TF_polII_integration", "03_SM_cluster_level_analysis")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_polIIDegs <- here::here("analysis", "06_polII_diff", "polII_DEGs.fpkm_filtered.combined.tab")
file_smRegions <- here::here("data", "reference_data", "SM_genes_regions.tab")
file_smRegionOrder <- here::here("data", "reference_data", "SM_regions_order.tab")
file_smBackbone <- here::here("data", "reference_data", "backbone_enz_type.txt")
file_smtfGridLayout <- paste(outDir, "/SMTF_grid_layout.tab", sep = "")

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
  dplyr::arrange(SM_ID, SMTF)

productionDataList <- purrr::transpose(productionData) %>% 
  purrr::set_names(nm = purrr::map(., "SMTF"))

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

combinedDegs <- suppressMessages(readr::read_tsv(file = file_polIIDegs))

## extract all SM cluster genes
smGeneInfo <- AnnotationDbi::select(
  x = orgDb, keys = keys(x = orgDb, keytype = "SM_ID"),
  columns = c("GID", "SM_CLUSTER"), keytype = "SM_ID"
) %>% 
  dplyr::rename(geneId = GID)


geneInfo <- AnnotationDbi::select(
  x = orgDb, keys = unique(smGeneInfo$geneId),
  columns = c("GENE_NAME", "DESCRIPTION"), keytype = "GID"
)


genePos <- as.data.frame(
  GenomicFeatures::genes(x = txDb, filter = list(gene_id = unique(smGeneInfo$geneId)))
) %>% 
  dplyr::left_join(
    y = geneInfo, by = c("gene_id" = "GID")
  ) %>% 
  dplyr::select(
    geneId = gene_id, geneName = GENE_NAME, chr = seqnames, start, end, strand
  )

smGeneInfo <- dplyr::left_join(x = smGeneInfo, y = genePos, by = "geneId") %>% 
  dplyr::arrange(chr, start, SM_ID)

clusterInfoList <- dplyr::select(smGeneInfo, SM_ID, SM_CLUSTER, geneId) %>% 
  split(.$SM_ID) %>% 
  purrr::map(.f = ~map(.x = ., .f = unique))


regionOrd <- suppressMessages(readr::read_tsv(file_smRegionOrder))

smGeneRegions <- suppressMessages(readr::read_tsv(file = file_smRegions)) %>% 
  dplyr::left_join(y = genePos, by = "geneId") %>%
  dplyr::arrange(chr, start, SM_IDs) %>% 
  dplyr::left_join(
    y = regionOrd, by = "regionId"
  )

geneTypes <- suppressMessages(readr::read_tsv(file = file_smBackbone)) %>% 
  dplyr::rename(geneSubtype = backbone_type) %>% 
  dplyr::mutate(geneType = "backbone")

geneTypes <- dplyr::select(productionData, geneId = SMTF) %>% 
  dplyr::mutate(geneType = "SMTF", geneSubtype = "SMTF") %>% 
  dplyr::bind_rows(geneTypes)

gridLayout <- suppressMessages(readr::read_tsv(file = file_smtfGridLayout))

###########################################################################
## prepare basic coordinates data for plotting

# dplyr::group_by(smGeneRegions, regionId) %>%
#   dplyr::mutate(
#     nGenes = n()
#   ) %>%
#   readr::write_tsv(file = paste(outPrefix, ".clusters.tab", sep = ""))

## 310 x 220
baseData <- dplyr::arrange(.data = smGeneRegions, chr, start) %>% 
  dplyr::left_join(y = geneTypes, by = "geneId") %>% 
  dplyr::mutate_at(
    .vars = c("chr", "geneId", "regionId"), .funs = ~forcats::as_factor(.)
  ) %>% 
  dplyr::group_by(regionId) %>%
  dplyr::arrange(start, .by_group = TRUE) %>%
  dplyr::mutate(
    genePos = 1:n()
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(
    chr = forcats::fct_drop(chr),
    regionId = forcats::fct_relevel(regionId, regionOrd$regionId)
  ) %>% 
  tidyr::replace_na(replace = list(geneType = "other", geneSubtype = "other"))

## X and Y coordinates for rectangular arrangement of the genes
baseData <- dplyr::group_by(baseData, side) %>%
  dplyr::arrange(regionId, chr, start, .by_group = TRUE) %>%
  dplyr::mutate(
    x = dplyr::case_when(
      side == "bottom" ~ as.integer(1:n()+ 5),
      side == "right" ~ as.integer(rep(220+10, n())),
      side == "top" ~ as.integer(n():1 + 5),
      side == "left" ~ as.integer(rep(1, n())),
    ),
    y = dplyr::case_when(
      side == "bottom" ~ as.integer(rep(1, n())),
      side == "right" ~ as.integer(1:n()+5),
      side == "top" ~ as.integer(rep(310+10, n())),
      side == "left" ~ as.integer(n():1+5),
    )
  ) %>% 
  dplyr::ungroup()

## function to convert the linear arrangement of points to evenly distributed
## points on bezier curve
points_on_curve <- function(df){
  
  ## no change for the clusters with <6 genes
  if(nrow(df) < 6){
    return(
      dplyr::mutate(df, xArc = x, yArc = y)
    )
  }
  
  curveControlPoint <- 1
  
  ## curve direction should change based on its side
  curveSign <- dplyr::case_when(
    df$side[1] == "top" ~ c(median(df$x), df$y[1] + quantile(1:nrow(df), curveControlPoint)),
    df$side[1] == "bottom" ~ c(median(df$x), df$y[1] - quantile(1:nrow(df), curveControlPoint)),
    df$side[1] == "right" ~ c(df$x[1] + quantile(1:nrow(df), curveControlPoint), median(df$y)),
    df$side[1] == "left" ~ c(df$x[1] - quantile(1:nrow(df), curveControlPoint), median(df$y)),
  )
  
  ## conver liner arrangement of points to bezier curve arrangement
  p <- matrix(
    data = c(
      df$x[1], df$y[1],
      curveSign,
      tail(df$x, 1), tail(df$y, 1)
    ),
    nrow=3, ncol=2, byrow=TRUE
  )
  
  pob <- pointsOnBezier(p=p, n=nrow(df), method="evenly_spaced", print.progress = T)
  # bp <- bezier(t=seq(0, 1, length=100), p=p)
  # plot(bp, cex=0.5, asp=1)
  # points(pob$points, col="red")
  
  df <- dplyr::mutate(
    df,
    xArc = pob$points[,1],
    yArc = pob$points[,2],
  )
  
  return(df)
}


# df <- dplyr::filter(baseData, regionId == "region_54")
# 
# dfNew <- points_on_curve(df = df)
# ggplot(data = dfNew, mapping = aes(xArc, yArc)) +
#   geom_point()

## convert the linear arrangement of genes to curved so point size can be increased
baseDataCurves <- dplyr::group_by(baseData, regionId) %>% 
  dplyr::group_modify(
    .f = ~points_on_curve(df = .)
  ) %>% 
  dplyr::ungroup()


## SMTF coordinates on the layout for drawing SMTF-target links
oetfCoord <- dplyr::select(productionData, geneId = SMTF, SMTF_name) %>% 
  dplyr::left_join(
    y = baseDataCurves, by = c("geneId")
  ) %>% 
  dplyr::select(
    OESMTF = geneId, SMTF_name, xArc.smtf = xArc, yArc.smtf = yArc,
    x.smtf = x, y.smtf = y, side.smtf = side, region.smtf = regionId
  ) %>% 
  dplyr::left_join(
    y = gridLayout, by = "OESMTF"
  )


## these genes will be labled with the cluster names
smClusterLabel <- dplyr::group_by(smGeneInfo, SM_ID) %>% 
  dplyr::arrange(chr, start, .by_group = TRUE) %>% 
  dplyr::mutate(
    n = n(),
    position = 1:n(),
    SM_CLUSTER = stringr::str_replace(
      string = SM_CLUSTER, pattern = "\\s*cluster", replacement = ""
    ),
    SM_CLUSTER = stringr::str_wrap(string = SM_CLUSTER, width = 15)
  ) %>% 
  dplyr::slice(n = ceiling(n()/2)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(chr, start) %>% 
  dplyr::select(geneId, SM_CLUSTER) %>% 
  dplyr::left_join(
    y = baseDataCurves, by = c("geneId")
  )


## outer layout of genes as curves arranged on rectangle perimeter
pt_layout <- ggplot2::ggplot(
) +
  geom_point(
    data = baseDataCurves,
    mapping = aes(x = xArc, y = yArc, color = geneType),
    shape = 16
  ) +
  scale_color_manual(
    values = c("SMTF" = "red", "backbone" = "blue", "other" = "black")
  ) +
  geom_text_repel(
    data = dplyr::filter(smClusterLabel, side == "left"),
    mapping = aes(x = xArc, y = yArc, label = SM_CLUSTER),
    hjust = 1, direction = "y", angle = 0, nudge_x = -10, fontface = "italic", size = 6
  ) +
  geom_text_repel(
    data = dplyr::filter(smClusterLabel, side == "right"),
    mapping = aes(x = xArc, y = yArc, label = SM_CLUSTER),
    hjust = 0, direction = "y", angle = 0, nudge_x = 10, fontface = "italic", size = 6
  ) +
  geom_text_repel(
    data = dplyr::filter(smClusterLabel, side == "top"),
    mapping = aes(x = xArc, y = yArc, label = SM_CLUSTER),
    hjust = 0, direction = "x", angle = 90, nudge_y = 10, fontface = "italic", size = 6
  ) +
  geom_text_repel(
    data = dplyr::filter(smClusterLabel, side == "bottom"),
    mapping = aes(x = xArc, y = yArc, label = SM_CLUSTER),
    hjust = 1, direction = "x", angle = 90, nudge_y = -10, fontface = "italic", size = 6
  ) +
  scale_x_continuous(expand = expansion(add = 50)) +
  scale_y_continuous(expand = expansion(add = 50)) +
  theme_nothing() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA)
  )


pdf(file = paste(outPrefix, ".layout.pdf", sep = ""), width = 14, height = 18)
pt_layout
dev.off()

png(file = paste(outPrefix, ".layout.png", sep = ""), width = 5000, height = 6500, res = 350)
pt_layout
dev.off()


# pdf(file = paste(outPrefix, ".rect_layout.pdf", sep = ""), width = 10, height = 14)
# 
# ggplot2::ggplot(
# ) +
#   geom_point(
#     data = baseDataCurves,
#     mapping = aes(x = x, y = y, color = OETF),
#     shape = 16
#   ) +
#   scale_color_manual(
#     values = c("yes" = "red", "no" = "black")
#   ) +
#   theme_nothing()
# 
# dev.off()

##################################################################################
## get DEG+binding data for each SM cluster genes from all SMTFOE experiments
rowId <- 1
mergedData <- NULL


for (rowId in 1:nrow(productionData)) {
  
  tfSampleId <- productionData$tfId[rowId]
  degId <- productionData$degId[rowId]
  
  tfSelfClusterGenes <- clusterInfoList[[
    productionDataList[[productionData$SMTF[rowId]]]$SM_ID
  ]]$geneId
  
  
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
  diffData <- dplyr::filter(combinedDegs, comparison == degId) %>%
    dplyr::select(geneId, comparison, !!col_lfc, shrinkLog2FC, pvalue, padj, maxFpkm, fpkmFilter)
  
  bindingDegData <- dplyr::left_join(
    x = baseDataCurves, y = diffData, by = "geneId"
  ) %>%
    dplyr::left_join(y = peakAn, by = "geneId") %>%
    dplyr::mutate(
      OESMTF = !!productionData$SMTF[rowId],
      OESMTF_name = !!productionData$SMTF_name[rowId]
    ) %>%
    dplyr::select(
      OESMTF, OESMTF_name, geneId, geneName, everything()
    )
  
  ## add self cluster target column
  if (!is.null(tfSelfClusterGenes)) {
    bindingDegData <- dplyr::left_join(
      x = bindingDegData,
      y = tibble::tibble(geneId = tfSelfClusterGenes, selfCluster = "yes"),
      by = "geneId"
    ) %>% 
      tidyr::replace_na(replace = list(selfCluster = "no"))
  } else{
    bindingDegData <- dplyr::mutate(
      bindingDegData,
      selfCluster = "no"
    )
  }
  
  mergedData <- dplyr::bind_rows(mergedData, bindingDegData)
  
}


## add filters and additional annotations  for the genes
mergedData2 <- dplyr::mutate(
  mergedData,
  significance = dplyr::case_when(
    !!sym(col_pval) <= cutoff_fdr & abs(!!sym(col_lfc)) >= cutoff_lfc & fpkmFilter == "pass" ~ "significant",
    TRUE ~ "non-significant"
  ),
  binding = dplyr::if_else(
    condition = hasPeak, true = "bound", false = "not-bound", missing = "not-bound"
  )
) %>% 
  dplyr::left_join(
    y = oetfCoord, by = c("OESMTF")
  ) %>% 
  dplyr::mutate(
    regulationType = dplyr::case_when(
      binding == "bound" & log2FoldChange < 0 & significance == "significant" ~ "bound+down",
      binding == "bound" & log2FoldChange > 0 & significance == "significant" ~ "bound+up",
      binding == "bound" ~ "bound",
      log2FoldChange < 0 & significance == "significant" ~ "down",
      log2FoldChange > 0 & significance == "significant" ~ "up",
      TRUE ~ "no-regulation"
    )
  ) %>% 
  dplyr::mutate(
    regulationType = forcats::fct_relevel(
      .f = regulationType,
      "bound+up", "bound+down", "up", "down", "bound", "no-regulation"
    ),
    binding = forcats::fct_relevel(
      .f = binding, "bound", "not-bound"
    ),
    selfCluster = forcats::fct_relevel(
      .f = selfCluster, "yes", "no"
    )
  ) %>% 
  dplyr::mutate(
    slope = (y - y.smtf)/(x - x.smtf),
    ## decide the curvature for geom_curve based on the slope of the link and direction
    curvature = dplyr::case_when(
      abs(slope) == Inf & side.smtf == "right" ~ -0.5*sign(slope),
      slope == 0 & side.smtf == "bottom" & x.smtf < x ~ -0.5,
      slope == 0 & side.smtf == "top" & x.smtf > x ~ -0.5,
      (side.smtf == "right" | side.smtf == "left") & slope < 0 ~ -0.5,
      (side.smtf == "top" | side.smtf == "bottom") & slope > 0 ~ -0.5,
      TRUE ~ 0.5
    ),
    curvature = dplyr::case_when(
      region.smtf == regionId ~ -1*curvature,
      TRUE ~ curvature
    )
  )


# table(mergedData2$significance, mergedData2$binding)
## AN7061 region_19 cluster_17  top
## AN7921 region_07 cluster_07  bottom
## AN7820	region_21	cluster_20	left
## AN0153	region_56	cluster_12	right

tmpMergedDf <- dplyr::filter(
  mergedData2,
  # OESMTF == "AN0153"
  OESMTF %in% c("AN7061", "AN7921", "AN7820", "AN0153")
)

# readr::write_tsv(x = tmpMergedDf, file = paste(outPrefix, ".test_data.tab", sep = ""))

## plotting TFOE regulation rectangle facets showing binding and DEG using links
pt_box_factets <- ggplot2::ggplot(
) +
  geom_segment(
    data = dplyr::filter(mergedData2, OESMTF != geneId) %>% 
      dplyr::filter(significance != "non-significant" | binding == "bound") %>% 
      dplyr::arrange(desc(regulationType)),
    mapping = aes(x = 230/2, y = 320/2, xend = xArc, yend = yArc, color = regulationType)
  ) +
  geom_point(
    data = dplyr::arrange(mergedData2, desc(selfCluster)),
    mapping = aes(x = xArc, y = yArc, color = selfCluster),
    shape = 16
  ) +
  geom_text(
    data = dplyr::filter(mergedData2, x == 115, side == "bottom"),
    mapping = aes(x = 115, y = yArc - 10, label = OESMTF_name),
    vjust = 1, size = 7, fontface = "bold.italic"
  ) +
  scale_color_manual(
    values = c("yes" = "green", "no" = "black",
               "bound" = "grey",
               "bound+down" = "blue", "down" = "#74add1",
               "bound+up" = "red", "up" = "#fdae61")
  ) +
  scale_x_continuous(expand = expansion(add = 5)) +
  scale_y_continuous(expand = expansion(add = c(50, 5))) +
  facet_grid(
    ## use facet grid for customized arrangement of panels
    rows = vars(gridRow), cols = vars(gridCol)
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.text = element_blank()
    # strip.text = element_text(size = 16, face = "bold", margin = margin(0,0,0,0))
  )


pdf(file = paste(outPrefix, ".SMTF_links_box.pdf", sep = ""), width = 10, height = 14)
pt_box_factets
dev.off()



pt_box_curves <- ggplot2::ggplot(
) +
  geom_curve(
    data = dplyr::filter(tmpMergedDf, OESMTF != geneId, curvature == 0.5) %>% 
      dplyr::filter(significance != "non-significant" | binding == "bound") %>%
      dplyr::arrange(desc(regulationType)),
    curvature = 0.5,
    mapping = aes(x = xArc.smtf, y = yArc.smtf, xend = xArc, yend = yArc, color = regulationType)
  ) +
  geom_curve(
    data = dplyr::filter(tmpMergedDf, OESMTF != geneId, curvature == -0.5) %>% 
      dplyr::filter(significance != "non-significant" | binding == "bound") %>%
      dplyr::arrange(desc(regulationType)),
    mapping = aes(x = xArc.smtf, y = yArc.smtf, xend = xArc, yend = yArc, color = regulationType),
    curvature = -0.5
  ) +  geom_point(
    data = dplyr::arrange(tmpMergedDf, desc(selfCluster)),
    mapping = aes(x = xArc, y = yArc, color = selfCluster),
    shape = 16
  ) +
  scale_color_manual(
    values = c("yes" = "green", "no" = "black",
               "bound" = "grey",
               "bound+down" = "blue", "down" = "#74add1",
               "bound+up" = "red", "up" = "#fdae61")
  ) +
  scale_x_continuous(expand = expansion(add = 5)) +
  scale_y_continuous(expand = expansion(add = 5)) +
  facet_wrap(
    facets = vars(OESMTF_name),
    strip.position = "bottom"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    # strip.text = element_blank()
    strip.text = element_text(size = 16, face = "bold", margin = margin(0,0,0,0))
  )


pt_box_curves





