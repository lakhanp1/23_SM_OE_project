suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(biovizBase))

## this script generate karyogram to show A. nidulans SM cluster genes on chromosomes

rm(list = ls())

##################################################################################
analysisName <- "SM_cluster_karyogram"
outDir <- here::here("analysis", "06_SM_cluster_summary")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
file_smBed <- "E:/Chris_UM/Database/A_Nidulans/SM_genes.bed"

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF


##################################################################################
productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::filter(has_polII_ChIP == "has_data", has_TF_ChIP == "has_data", copyNumber == "sCopy")



smInfo <- suppressMessages(
  AnnotationDbi::select(
    x = orgDb, keys = keys(orgDb, keytype = "SM_ID"),
    columns = c("GID", "SM_ID", "SM_CLUSTER"), keytype = "SM_ID")
) %>% 
  dplyr::rename(geneId = GID)

bgcWithData <- unique(smInfo$SM_ID[which(smInfo$geneId %in% productionData$geneId)])

smInfo$hasData <- smInfo$SM_ID %in% bgcWithData

smGr <- GenomicFeatures::genes(
  x = txDb, filter = list(gene_id = smInfo$geneId)
) %>% 
  as.data.frame %>% 
  dplyr::left_join(y = smInfo, by = c("gene_id" = "geneId")) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

seqlevels(smGr) <- seqlevels(txDb)
seqinfo(smGr) <- seqinfo(txDb)

# smGr <- keepSeqlevels(x = smGr, value = sort(seqlevels(smGr)))

pt_theme <- theme(
  plot.title = element_text(size = 18, face = "bold", hjust = 1),
  axis.text = element_text(size = 16, face = "bold", color = "black"),
  axis.ticks.length = unit(2, "mm"),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title = element_blank(),
  legend.position = "bottom",
  legend.text = element_text(size = 16, face = "bold"),
  strip.background = element_rect(fill = "white"),
  strip.text.y = element_text(angle = 0),
  panel.border = element_blank(),
  panel.grid = element_blank()
)

pt_all <- autoplot(object = seqinfo(smGr), layout = "karyogram") +
  ggbio::layout_karyogram(
    data = smGr, geom = "rect", color = "black", 
    mapping = aes(fill = "black"),
    ylim = c(2, 8)
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("black" = "black"),
    labels = c("black" = "BGCs")
  ) +
  labs(title = "Distribution of secondary metabolism gene clusters in A. nidulans genome") +
  theme_clear() +
  pt_theme

pdf(file = paste(outPrefix, ".pdf", sep = ""), width = 12, height = 6)
pt_all
dev.off()



pt_data <- autoplot(object = seqinfo(smGr), layout = "karyogram") +
  ggbio::layout_karyogram(
    data = smGr, geom = "rect", 
    mapping = aes(color = hasData, fill = hasData),
    ylim = c(2, 8)
  ) +
  scale_color_manual(
    name = NULL,
    values = c("TRUE" = "blue", "FALSE" = "black"),
    labels = c("TRUE" = "has data", "FALSE" = "no data")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("TRUE" = "blue", "FALSE" = "black"),
    labels = c("TRUE" = "has data", "FALSE" = "no data")
  ) +
  labs(title = "Distribution of secondary metabolism gene clusters in A. nidulans genome") +
  theme_clear() +
  pt_theme


pdf(file = paste(outPrefix, ".data.pdf", sep = ""), width = 12, height = 6)
pt_data
dev.off()

