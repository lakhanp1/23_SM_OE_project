library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(TxDb.Anidulans.FGSCA4.AspGD.GFF)
library(here)
library(ChIPpeakAnno)
library(ChIPseeker)


## 1) annotate peaks
## 2) create gene level peak annotation data

rm(list = ls())

##################################################################################

file_exptInfo <- here::here("data", "referenceData/sample_info.txt")

file_genes <- here::here("data", "referenceData/AN_genes_for_polII.bed")
orgDb <- org.Anidulans.FGSCA4.eg.db
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

##################################################################################
## ChIPseeker test
peakAn <- annotatePeak(
  peak = peaksGr,
  tssRegion = c(-1000, 10),
  TxDb = txDb,
  ignoreDownstream = TRUE
)

plotAnnoBar(peakAn)
upsetplot(peakAn)

peakAnDf <- as.data.frame(peakAn)

readr::write_tsv(x = peakAnDf, path = "peak_annotation.chipseeker.tab")

##################################################################################
## ChIPpeakAnno test

annoData <- toGRanges(data = txDb, feature = "gene")

annotatedPeak <- annotatePeakInBatch(
  myPeakList = peaksGr,
  AnnotationData = annoData,
  output = "both",
  bindingRegion = c(-1000, 100),
  select = "all"
  )

readr::write_tsv(x = as.data.frame(annotatedPeak), path = "peak_annotation.chippeakanno.tab")

bdpAnn <- annotatePeakInBatch(
  myPeakList = peaksGr,
  AnnotationData = annoData,
  output = "nearestBiDirectionalPromoters",
  bindingRegion = c(-1000, 100),
  select = "all"
)



