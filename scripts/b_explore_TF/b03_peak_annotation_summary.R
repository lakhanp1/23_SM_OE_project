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

analysisName <- "peak_annotation"
outDir <- here::here("analysis", "09_TF_binding", "01_peak_annotation")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_dataSummary <- here::here("data", "reference_data", "production_data.summary.tab")
file_tfSamples <- here::here("data", "reference_data", "production_data.tf_samples.txt")
file_peakSummary <- here::here("analysis", "02_QC_TF", "TF_ChIP_summary.best_replicates.tab")
TF_dataPath <- here::here("data", "TF_data")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

cutoff_macs2Pval <- 20

##################################################################################

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

# dataSummary <- suppressMessages(readr::read_tsv(file = file_dataSummary)) %>% 
#   dplyr::filter(has_TF_ChIP == "has_data", has_polII_ChIP == "has_data")

tfSampleList <- suppressMessages(readr::read_tsv(file = file_tfSamples))


tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSampleList$sampleId,
  dataPath = TF_dataPath)

# glimpse(tfInfo)

clusterInfo <- AnnotationDbi::select(
  x = orgDb, keys = tfInfo$SM_TF, columns = c("GENE_NAME"), keytype = "GID"
) %>% 
  dplyr::mutate(
    geneLabel = paste(GID, " (", GENE_NAME, ")", sep = ""),
    geneLabel = if_else(condition = GID == GENE_NAME, true = GID, false = geneLabel)
  )

tfInfo <- dplyr::left_join(x = tfInfo, y = clusterInfo, by = c("SM_TF" = "GID"))

tfInfoList <- purrr::transpose(tfInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

tfSummary <- suppressMessages(readr::read_tsv(file = file_peakSummary))

genesDf <- as.data.frame(GenomicFeatures::genes(x = txDb)) %>% 
  dplyr::select(geneId = gene_id, strand)

##################################################################################

rowId <- 1

mergedAn <- NULL

for (rowId in 1:nrow(tfInfo)) {
  sampleId <- tfInfo$sampleId[rowId]
  
  anDf <- suppressMessages(readr::read_tsv(tfInfoList[[sampleId]]$peakAnno)) %>% 
    dplyr::mutate(
      sampleId = sampleId,
      SM_TF = tfInfo$SM_TF[rowId],
      geneLabel = tfInfo$geneLabel[rowId]
    ) %>% 
    dplyr::filter(peakPval >= cutoff_macs2Pval) %>% 
    dplyr::left_join(y = genesDf, by = "geneId")
  
  mergedAn <- dplyr::bind_rows(mergedAn, anDf)
  
}


outDf <- dplyr::filter(mergedAn, peakCategory != "intergenic" | peakCategory != "blacklist") %>% 
  dplyr::select(
    sampleId, SM_TF, peakId, peakEnrichment, peakPval, peakRegion, geneId, peakDist,
    peakAnnotation, peakCategory, bidirectional, relativeSummitPos, relativePeakPos
  )

readr::write_tsv(x = outDf, file = paste(outPrefix, ".combined.tab", sep = ""))

table(mergedAn$peakAnnotation, mergedAn$peakCategory)



combinedAn <- dplyr::mutate(
  mergedAn,
  newCategory = dplyr::case_when(
    peakAnnotation == "5UTR" ~ "promoter",
    peakAnnotation == "tx_start" ~ "promoter",
    peakAnnotation == "3UTR" ~ "nearEnd",
    peakAnnotation == "tx_end" ~ "nearEnd",
    peakAnnotation == "EXON" ~ "Exon/Intron",
    peakAnnotation == "INTRON" ~ "Exon/Intron",
    peakAnnotation == "include_tx" ~ "featureInPeak",
    peakAnnotation == "upstream_intergenic" ~ "intergenic",
    TRUE ~ peakAnnotation
  )
)


combinedAn <- dplyr::mutate(
  combinedAn,
  peakAnnotation = forcats::fct_relevel(
    peakAnnotation, "upstream_intergenic", "upstream", "promoter", "5UTR", "tx_start",
    "EXON", "INTRON", "include_tx", "3UTR", "tx_end", "intergenic"
  ),
  peakCategory = forcats::fct_relevel(
    peakCategory, "upstreamTss", "nearStart", "peakInFeature", 
    "featureInPeak", "nearEnd", "intergenic"
  ),
  newCategory = forcats::fct_relevel(
    newCategory, "upstream", "promoter",
    "Exon/Intron", "featureInPeak", "nearEnd", "intergenic"
  )
)


peakTypeColors <- structure(
  RColorBrewer::brewer.pal(n = 7, name = "Paired"),
  names = levels(combinedAn$newCategory))


pt_anBar <- ggplot(data = combinedAn) +
  geom_bar(mapping = aes(y = geneLabel, fill = newCategory),
           position = position_fill(reverse = TRUE)) +
  scale_x_continuous(
    labels = scales::percent, expand = expansion(add = 0)
  ) +
  # guides(fill = guide_legend(title.position = "bottom", title.hjust = 0.5)) +
  scale_fill_manual(name = "Peak annotation", values = peakTypeColors) +
  labs(
    title = "SM TFOE ChIPseq peak annotation distribution"
  ) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_blank(),
    title = element_text(size = 18),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 16, face = "bold"),
    plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
    panel.grid = element_blank()
  )

png(filename = paste(outPrefix, ".distribution.png", sep = ""), width = 3000, height = 3000, res = 300)
pt_anBar
dev.off()









