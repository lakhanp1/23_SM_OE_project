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
outDir <- here::here("analysis", "09_TF_binding", "peak_annotation")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_dataSummary <- here::here("data", "reference_data", "raw_data_summary.tab")
# file_tfSamples <- here::here("data", "reference_data", "tf_sample_ids.best_rep.txt")
file_peakSummary <- here::here("analysis", "02_QC_TF", "TF_ChIP_summary.best_replicates.tab")
TF_dataPath <- here::here("data", "TF_data")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

cutoff_macs2Pval <- 20

##################################################################################

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

dataSummary <- suppressMessages(readr::read_tsv(file = file_dataSummary)) %>% 
  dplyr::filter(has_TF_ChIP == "has_data", has_polII_ChIP == "has_data")


tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = dataSummary$tfId,
  dataPath = TF_dataPath)

# glimpse(tfInfo)

tfInfoList <- purrr::transpose(tfInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

tfSummary <- suppressMessages(readr::read_tsv(file = file_peakSummary))

genesDf <- as.data.frame(GenomicFeatures::genes(x = txDb)) %>% 
  dplyr::select(geneId = gene_id, strand)

##################################################################################

i <- 1

mergedAn <- NULL

for (i in 1:nrow(tfInfo)) {
  sampleId <- tfInfo$sampleId[i]
  
  anDf <- suppressMessages(readr::read_tsv(tfInfoList[[sampleId]]$peakAnno)) %>% 
    dplyr::mutate(sampleId = sampleId) %>% 
    dplyr::filter(peakPval >= cutoff_macs2Pval) %>% 
    dplyr::left_join(y = genesDf, by = "geneId") %>% 
    dplyr::left_join(
      y = dplyr::select(dataSummary, SM_TF = geneId, tfId), by = c("sampleId" = "tfId")
    )
  
  mergedAn <- dplyr::bind_rows(mergedAn, anDf)
  
}


outDf <- dplyr::filter(mergedAn, peakCategory != "intergenic") %>% 
  dplyr::select(
    sampleId, SM_TF, peakId, peakEnrichment, peakPval, peakRegion, geneId, peakDist,
    peakAnnotation, peakCategory, bidirectional, relativeSummitPos, relativePeakPos
  )

readr::write_tsv(x = outDf, path = paste(outPrefix, ".combined.tab", sep = ""))

table(mergedAn$peakAnnotation, mergedAn$peakCategory)



combinedAn <- dplyr::mutate(
  mergedAn,
  newCategory = dplyr::case_when(
    peakAnnotation == "5UTR" ~ "nearStart",
    peakAnnotation == "tx_start" ~ "nearStart",
    peakAnnotation == "3UTR" ~ "nearEnd",
    peakAnnotation == "tx_end" ~ "nearEnd",
    peakAnnotation == "EXON" ~ "Exon/Intron",
    peakAnnotation == "INTRON" ~ "Exon/Intron",
    peakAnnotation == "include_tx" ~ "featureInPeak",
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
    newCategory, "upstream_intergenic", "upstream", "promoter", "nearStart",
    "Exon/Intron", "featureInPeak", "nearEnd", "intergenic"
  )
)


peakTypeColors <- structure(
  RColorBrewer::brewer.pal(n = 8, name = "Paired"),
  names = levels(combinedAn$newCategory))


pt_anBar <- ggplot(data = combinedAn) +
  geom_bar(mapping = aes(y = SM_TF, fill = newCategory),
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









