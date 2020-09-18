suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(BSgenome.Anidulans.FGSCA4.AspGD))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggpubr))


## integrate polII DEG data and motif search data for AN0153

rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

##################################################################################

analysisName <- "AN0153_DEG_motifs_hits"
outDir <- here("analysis", "04_AN0153_analysis", "01_motif_enrichment")

file_fimo <- here("analysis", "04_AN0153_analysis", "01_motif_enrichment",
                  "fimo.AN0153_OE_vs_WT_DEG", "fimo.tsv")

outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")

TF_dataPath <- here::here("data", "TF_data")
diffDataPath <- here::here("analysis", "06_polII_diff")

samples_tf <- "AN0153_sCopy_OE_16h_HA_ChIPMix64_3"
polIIDiffIds <- "AN0153_sCopy_OE_vs_MH11036"


if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

orgDb <- org.Anidulans.FGSCA4.eg.db
genome <- BSgenome.Anidulans.FGSCA4.AspGD
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF


cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- -1 * cutoff_lfc

col_lfc <- "log2FoldChange"
col_fdr <- "padj"

##################################################################################

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison == polIIDiffIds)

exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = samples_tf,
                                   dataPath = TF_dataPath)

peakAn <- import_peak_annotation(
  sampleId = exptData$sampleId,
  peakAnnoFile = exptData$peakAnno,
  renameColumn = FALSE
  ) %>% 
  dplyr::filter(peakPosition == "TSS" & peakDist >= -1000 & peakPval >= 20) %>% 
  dplyr::group_by(geneId) %>% 
  dplyr::arrange(desc(peakDist)) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(peakId, peakEnrichment, peakPval, peakDist, geneId) %>% 
  dplyr::mutate(hasPeak = "bound")



diffData <- suppressMessages(
  readr::read_tsv(file = rnaseqInfo$deg)
)

significantDegs <-  dplyr::filter(
  diffData,
  abs(!!sym(col_lfc)) >= cutoff_lfc & !!sym(col_fdr) <= cutoff_fdr
)

fimo <- suppressMessages(readr::read_tsv(file = file_fimo, comment = "#"))

geneList <- split(x = significantDegs$geneId,
                  f = forcats::as_factor(significantDegs$diff_l2fc))

fimoSummary <- dplyr::group_by(fimo, sequence_name) %>% 
  dplyr::summarise(nMotifs = n()) %>% 
  dplyr::rename(geneId = sequence_name)

degGr <- GenomicFeatures::genes(x = txDb, filter = list(gene_id = significantDegs$geneId))

degTss <- GenomicRanges::resize(x = degGr, width = 1, fix = "start")

degTss$name <- degTss$gene_id

plotDegs <- dplyr::left_join(
  x = tibble(geneId = degTss$gene_id),
  y = dplyr::select(significantDegs, geneId, log2FoldChange, diff_l2fc),
  by = "geneId"
) %>% 
  dplyr::left_join(y = fimoSummary, by = "geneId") %>% 
  dplyr::left_join(y = peakAn, by = "geneId") %>% 
  dplyr::mutate(diff_l2fc = forcats::fct_relevel(.f = diff_l2fc, "up")) %>% 
  tidyr::replace_na(
    list(nMotifs = 0, peakEnrichment = 0, peakPval = 0, hasPeak = "notBound")
  ) %>% 
  tidyr::unite(col = "cluster", diff_l2fc, hasPeak, sep = ":", remove = FALSE) %>% 
  dplyr::mutate(
    cluster = forcats::fct_relevel(
      .f = cluster, "up:bound", "up:notBound", "down:bound", "down:notBound")
  )


geneStats <- table(plotDegs$cluster)
map <- structure(.Data = names(geneStats),
                 names =paste(names(geneStats), " (", geneStats, ")", sep = ""))

plotDegs$cluster <- forcats::fct_recode(.f = plotDegs$cluster, !!!map)

dplyr::glimpse(plotDegs)

readr::write_tsv(x = plotDegs, path = paste(outPrefix, ".heatmap.data.tab", sep = ""))

## plotting
mat <- chipmine::bigwig_profile_matrix(
  bwFile = exptData$bwFile,
  regions = degTss,
  extend = c(1000, 200),
  signalName = exptData$sampleId,
  targetName = "TSS")


quantile(mat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
tfColor <- colorRamp2(
  breaks = quantile(mat, c(0.50, 0.98), na.rm = T),
  colors =  c("black", "yellow")
)

anBar <- HeatmapAnnotation(
  motif_count = anno_barplot(
    x = plotDegs$nMotifs, width = unit(3, "cm")),
  width = unit(3, "cm"),
  which = "row",
  annotation_name_side = "top",
  annotation_label = "Motif",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 14),
  show_legend = TRUE
)

pt1 <- chipmine::profile_heatmap(
  profileMat = mat,
  signalName = "AN0153 binding",
  profileColor = tfColor,
  geneGroups = dplyr::select(plotDegs, geneId, cluster),
  showAnnotation = FALSE,
  column_title_gp = gpar(fontsize = 14),
  width = unit(10, "cm"),
  posLineGpar = gpar(col = "white", lty = 2, alpha = 0.4, lwd = 1),
  row_title_gp = gpar(fontsize = 18),
  heatmap_legend_param = list(
    legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
    labels_gp = gpar(fontsize = 16),
    title_gp = gpar(fontsize = 16, fontface = "bold")
  )
)

ht_lfc <- Heatmap(
  matrix = matrix(data = plotDegs$log2FoldChange, dimnames = list(plotDegs$geneId, "polII_lfc")),
  col = colorRamp2(breaks = c(-2, -0.5, 0.5, 2), colors = c("#313695", "white", "white", "#d73027")),
  name = "polII_lfc",
  column_title = "OE/WT DEG",
  show_column_names = FALSE,
  cluster_rows = FALSE, show_row_names = FALSE,
  width = unit(3, "cm"),
  heatmap_legend_param = list(
    legend_height = unit(4, "cm"), grid_width = unit(1, "cm"),
    labels_gp = gpar(fontsize = 16),
    title_gp = gpar(fontsize = 16, fontface = "bold")
  )
)

htList <- pt1$heatmap + ht_lfc + anBar

pdf(file = paste(outPrefix, ".heatmap.pdf", sep = ""), width = 12, height = 12)

draw(
  object = htList,
  column_title = "AN0153 polII DEGs, binding and motif counts",
  column_title_gp = gpar(fontsize = 20, fontface = "bold")
)

dev.off()


ggTheme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold"),
  )

gg_nMotif <- ggplot(data = plotDegs) +
  geom_boxplot(mapping = aes(x = cluster, y = nMotifs),
               size = 1, fill='#A4A4A4') +
  labs(title = "number of Tbf1 motifs in each DEG promoter (-500:TSS:100)") +
  coord_flip() +
  ggTheme

gg_peakEnrich <- ggplot(data = plotDegs) +
  geom_boxplot(mapping = aes(x = cluster, y = peakEnrichment),
               size = 1, fill='#E69F00') +
  labs(title = "macs2 peak enrichment for Tbf1 peaks") +
  coord_flip() +
  ggTheme

pdf(file = paste(outPrefix, ".boxplots.pdf", sep = ""), width = 14, height = 8)
ggpubr::ggarrange(gg_nMotif, gg_peakEnrich, ncol = 2, align = "h")
dev.off()




