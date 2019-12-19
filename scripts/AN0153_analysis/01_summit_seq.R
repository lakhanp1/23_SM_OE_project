library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(TxDb.Anidulans.AspGD.GFF)
library(BSgenome.Anidulans.AspGD.FGSCA4)
library(here)


## get summit sequence

rm(list = ls())

##################################################################################

analysisName <- "AN0153_summit_seq"
outDir <- here("analysis", "04_AN0153_analysis", "motif_enrichment")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

orgDb <- org.Anidulans.FGSCA4.eg.db
genome <- BSgenome.Anidulans.AspGD.FGSCA4

file_exptInfo <- here::here("data", "referenceData/sample_info.txt")

file_sampleIds <- here("analysis", "04_AN0153_analysis", "sample_ids.txt")

TF_dataPath <- here::here("data", "TF_data")

##################################################################################

sampleIds <- suppressMessages(readr::read_tsv(file = file_sampleIds, comment = "#"))

exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = sampleIds$sampleId,
                                   dataPath = TF_dataPath)

i <- 1

summit200 <- get_peak_summit_seq(file = exptData$peakFile[i],
                                 peakFormat = "narrowPeak",
                                 sampleId = exptData$sampleId[i],
                                 genome = genome, length = 200,
                                 column_name_prefix = FALSE)

summit200 <- dplyr::select(summit200, name, summitSeq) %>% 
  dplyr::rename(summitSeq200 = summitSeq)

summit500 <- get_peak_summit_seq(file = exptData$peakFile[i],
                                 peakFormat = "narrowPeak",
                                 sampleId = exptData$sampleId[i],
                                 genome = genome, length = 500,
                                 column_name_prefix = FALSE)

summit500 <- dplyr::select(summit500, name, summitSeq) %>% 
  dplyr::rename(summitSeq500 = summitSeq)


peakDf <- import_peaks_as_df(file = exptData$peakFile[i],
                             sampleId = exptData$sampleId[i],
                             peakFormat = "narrowPeak", rename = F)

peakDf <- dplyr::left_join(x = peakDf, y = summit200, by = c("peakId" = "name")) %>% 
  dplyr::left_join(y = summit500, by = c("peakId" = "name"))

readr::write_tsv(x = peakDf,
                 path = paste(outDir, "/", exptData$sampleId[i], ".summitSeq.tab", sep = ""))















