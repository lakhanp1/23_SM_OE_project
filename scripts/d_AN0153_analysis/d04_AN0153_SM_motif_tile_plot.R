suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(BSgenome.Anidulans.FGSCA4.AspGD))
suppressPackageStartupMessages(library(here))


## tile plot to visualize AN0153 binding, polII fold change and TTAGGG motif

rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

##################################################################################

analysisName <- "tfb1_binding.polII.motif"
outDir <- here("analysis", "04_AN0153_analysis", "01_motif_enrichment")

file_motifMatch <- here("analysis", "04_AN0153_analysis", "01_motif_enrichment",
                        "TTAGG_SM_gene_match.tab")

outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")

TF_dataPath <- here::here("data", "TF_data")
diffDataPath <- here::here("analysis", "07_polII_diff")

samples_tf <- "AN0153_sCopy_OE_16h_HA_ChIPMix46_1"
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

motifHits <- suppressMessages(readr::read_tsv(file = file_motifMatch))


diffData <- suppressMessages(
  readr::read_tsv(file = rnaseqInfo$deg)
)

smGenes <- AnnotationDbi::select(x = orgDb, keys = keys(x = orgDb, keytype = "SM_GENE"),
                                 columns = c("SM_CLUSTER", "SM_ID"), keytype = "SM_GENE") %>% 
  dplyr::rename(geneId = SM_GENE) %>% 
  dplyr::mutate(
    SM_ID = stringr::str_replace(string = SM_ID, pattern = "cluster_", replacement = ""),
    SM_CLUSTER = paste("(", SM_ID, ") ", SM_CLUSTER, sep = "")
  )

genesGrDf <- as.data.frame(GenomicFeatures::genes(x = txDb)) %>% 
  dplyr::select(seqnames, start, end, strand, geneId = gene_id)


tfBindingDf <- suppressMessages(readr::read_tsv(file = exptData$peakAnno)) %>% 
  dplyr::select(geneId, peakId, peakPval, peakDist, peakType, peakPosition, peakCategory) %>% 
  dplyr::filter(peakDist > -1000) %>%
  dplyr::group_by(geneId) %>% 
  dplyr::summarise(peakPosition = paste(sort(unique(peakPosition)), collapse = ","))

##################################################################################


plotDf <- dplyr::left_join(x = smGenes, y = motifHits, by = "geneId") %>% 
  dplyr::left_join(y = dplyr::select(diffData, geneId, diff_l2fc), by = "geneId") %>% 
  dplyr::left_join(y = genesGrDf, by = "geneId") %>% 
  dplyr::left_join(y = tfBindingDf, by = "geneId") %>% 
  dplyr::group_by(SM_ID) %>% 
  dplyr::arrange(start, .by_group = TRUE) %>% 
  dplyr::mutate(
    index = row_number(),
    midGeneIdx = round(n()/2)
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(
    centeredIndex = index + max(index)/2 - midGeneIdx,
    SM_CLUSTER = forcats::as_factor(SM_CLUSTER)
  )

glimpse(plotDf)

ptTheme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        # strip.text.x = element_text(hjust = 0.5, size = 12, face = "bold"),
        # strip.background = element_rect(fill="white", size = 0.2),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 13, face = "bold"),
        plot.margin = unit(rep(0.5, 4), "cm"))

xCol <- "centeredIndex"

pt <- ggplot() +
  geom_point(
    data = dplyr::filter(plotDf, !is.na(TTAGGG)),
    mapping = aes(x = !!sym(xCol), y = "motif", color = TTAGGG),
    size = 2, shape = 16
  ) +
  geom_tile(
    data = plotDf,
    mapping = aes(x = !!sym(xCol), y = "l2fc", fill = diff_l2fc),
    color = "black", size = 0.5, height = 1
  ) +
  geom_point(
    data = dplyr::filter(plotDf, !is.na(peakPosition)),
    mapping = aes(x = !!sym(xCol), y = "peak"), color = "black",
    size = 2, shape = 17
  ) +
  geom_text(
    ## adding geom_text() to move strip.text inside the box as label
    data = dplyr::group_by(plotDf, SM_CLUSTER) %>% 
      dplyr::arrange(index) %>% 
      dplyr::slice(n()) %>% 
      dplyr::ungroup(),
    mapping = aes(x = max(plotDf$index)/2, y = "label", label = SM_CLUSTER),
    size = 5, fontface = "bold"
  ) +
  labs(
    title = "AN0153 OE peak, OE_polII/WT_polII and TTAGGG motif on SM clusters",
    subtitle = "row1: AN0153 ChIPseq peak || row2: polII AN0153:OE/WT DEG || row3: TTAGGG motif"
  ) +
  scale_fill_manual(
    values = c("down" = "blue", "noDEG" = "#F7F7F7", "up" = "red"),
  ) +
  scale_color_manual(
    values = c("promoter" = "#B35806", "3UTR" = "#542788"),
    breaks = c("promoter", "3UTR")
  ) +
  scale_x_continuous(expand = expansion(add = c(0.5, 0.5))) +
  scale_y_discrete(
    limits = c("motif", "l2fc", "peak", "label"),
    expand = expansion(add = c(0.5, 0.75))
  ) +
  facet_wrap(facets = SM_CLUSTER ~ ., scales = "free_y", ncol = 3, dir = "v") +
  ptTheme


pdf(file = paste(outPrefix, ".tile_plot.pdf", sep = ""), width = 12, height = 20)
pt
dev.off()











