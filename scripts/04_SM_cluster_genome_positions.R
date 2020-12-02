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

file_smBed <- "E:/Chris_UM/Database/A_Nidulans/SM_genes.bed"

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF



smGr <- rtracklayer::import(con = file_smBed, format = "bed")
seqlevels(smGr) <- seqlevels(txDb)

seqinfo(smGr) <- seqinfo(txDb)

smGr <- keepSeqlevels(x = smGr, value = sort(seqlevels(smGr)))

pt_theme <- theme(
  plot.title = element_text(size = 18, face = "bold"),
  axis.text = element_text(size = 16),
  axis.ticks.length = unit(2, "mm"),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title = element_blank(),
  strip.background = element_rect(fill = "white"),
  strip.text.y = element_text(angle = 0),
  panel.border = element_blank(),
  panel.grid = element_blank()
)

pt <- autoplot(object = seqinfo(smGr), layout = "karyogram") +
  layout_karyogram(data = smGr, geom = "rect", color = "black", fill = "black",
                   ylim = c(2, 8)) +
  labs(title = "Distribution of secondary metabolism gene clusters in A. nidulans genome") +
  # scale_x_continuous(
  #   expand = expansion(mult = 0.01)
  # ) +
  theme_clear() +
  pt_theme

pdf(file = paste(outPrefix, ".pdf", sep = ""), width = 12, height = 6)
pt
dev.off()







