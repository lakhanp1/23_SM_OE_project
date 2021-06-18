suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
# suppressPackageStartupMessages(library(TxDb.Anidulans.tRNA.removed))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ChIPpeakAnno))
suppressPackageStartupMessages(library(ChIPseeker))


## 1) annotate peaks using ChIPmine, ChIPseeker and ChIPpeakAnno
## 2) create gene level peak annotation data
## 3) compare the gene annotation

rm(list = ls())

##################################################################################

file_exptInfo <- here::here("data", "reference_data/sample_info.txt")

file_genes <- here::here("data", "reference_data/AN_genes_for_polII.bed")
orgDb <- org.Anidulans.FGSCA4.eg.db
# txDb <- TxDb.Anidulans.tRNA.removed
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")
other_dataPath <- here::here("data", "other_data")

file_tf_macs2 <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")
file_tf <- paste(TF_dataPath, "/", "sample_tf.list", sep = "")

##################################################################################

tfSampleList <- readr::read_tsv(file = file_tf_macs2, col_names = c("sampleId"),  comment = "#")

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSampleList$sampleId,
  dataPath = TF_dataPath)

i <- 3

peaksGr <- rtracklayer::import(con = tfInfo$peakFile[i], format = "narrowPeak")


txInfo <- suppressMessages(
  AnnotationDbi::select(
    x = txDb, keys = AnnotationDbi::keys(x = txDb, keytype = "TXID"),
    columns = c("GENEID", "TXNAME", "TXTYPE"), keytype = "TXID")) %>%
  dplyr::mutate(TXID = as.character(TXID)) %>%
  dplyr::rename(geneId = GENEID, txName = TXNAME, txType = TXTYPE)

txInfo <- dplyr::filter(txInfo, !txType %in% c("tRNA", "rRNA", "snRNA", "snoRNA")) %>% 
  dplyr::filter(!grepl(pattern = "uORF", x = geneId))

## extract transcript GRanges
transcriptsGr <- GenomicFeatures::transcripts(
  x = txDb, filter = list(tx_id = txInfo$TXID),
  columns = c("tx_id", "tx_name", "tx_type", "gene_id")
)

##################################################################################
## ChIPseeker test
peakAn <- annotatePeak(
  peak = peaksGr,
  tssRegion = c(-1000, 10),
  TxDb = txDb,
  # TxDb = transcriptsGr,
  ignoreDownstream = TRUE,
  overlap = "all",
  verbose = FALSE
)

plotAnnoBar(peakAn)
upsetplot(peakAn)

peakAnDf <- as.data.frame(peakAn)
# glimpse(peakAnDf)

readr::write_tsv(x = peakAnDf, path = "peak_annotation.chipseeker.tab")

##################################################################################
## ChIPpeakAnno test

# annoData <- toGRanges(data = txDb, feature = "gene")
# 
# annotatedPeak <- annotatePeakInBatch(
#   myPeakList = peaksGr,
#   AnnotationData = annoData,
#   output = "both",
#   bindingRegion = c(-1000, 100),
#   select = "all"
# )
# 
# annotatedPeak <- as.data.frame(annotatedPeak)
# 
# readr::write_tsv(x = annotatedPeak, path = "peak_annotation.chippeakanno.tab")
# 
# bdpAnn <- annotatePeakInBatch(
#   myPeakList = peaksGr,
#   AnnotationData = annoData,
#   output = "nearestBiDirectionalPromoters",
#   bindingRegion = c(-1000, 100),
#   select = "all"
# )

##################################################################################
## chipmine annotation
chipmineAn <- suppressMessages(readr::read_tsv(file = tfInfo$peakAnno[i]))

cmAn <- dplyr::select(chipmineAn, peakId, geneId, peakType, peakCategory, bidirectional,
                      peakPosition, relativePeakPos, peakDist, targetOverlap) %>% 
  dplyr::mutate(cm_an = TRUE)


# cpaAn <- dplyr::select(annotatedPeak, name, cpa_feature = feature,
#                        cpa_insideFeature = insideFeature,
#                        cpa_distancetoFeature = distancetoFeature,
#                        cpa_distancetoFeature = distancetoFeature) %>% 
#   dplyr::mutate(cpa_an = TRUE)

csAn <- dplyr::select(
  peakAnDf, name, cs_geneId = geneId,
  cs_distanceToTSS = distanceToTSS, cs_annotation = annotation
) %>% 
  dplyr::mutate(cs_an = TRUE)

compareDf <- dplyr::full_join(x = cmAn, y = csAn,
                              by = c("peakId" = "name", "geneId" = "cs_geneId"))


compareDf <- as.data.frame(peaksGr) %>% 
  dplyr::select(seqnames, start, end, peakId = name) %>% 
  dplyr::left_join(y = compareDf, by = "peakId") %>% 
  dplyr::mutate(
    setDiff = paste(cm_an, ":", cs_an, sep = ""),
    peakGene = paste(peakId, ":", geneId, sep = "")
  ) %>% 
  dplyr::select(
    seqnames, start, end, intersect(colnames(compareDf), colnames(cmAn)),
    setDiff, intersect(colnames(compareDf), colnames(csAn)), peakGene
  )

# glimpse(compareDf)

readr::write_tsv(x = compareDf, path = "peak_annotation.compare_cm_cs.tab")





