suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(BSgenome.Anidulans.FGSCA4.AspGD))
suppressPackageStartupMessages(library(here))


## get summit sequence for AN0153 ChIPseq peaks

rm(list = ls())

##################################################################################

analysisName <- "AN0153_summit_seq"
outDir <- here("analysis", "04_AN0153_analysis", "01_motif_enrichment")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

orgDb <- org.Anidulans.FGSCA4.eg.db
genome <- BSgenome.Anidulans.FGSCA4.AspGD

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")

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

summit200 <- dplyr::select(summit200, peakId, summitSeq) %>% 
  dplyr::rename(summitSeq200 = summitSeq)

summit500 <- get_peak_summit_seq(file = exptData$peakFile[i],
                                 peakFormat = "narrowPeak",
                                 sampleId = exptData$sampleId[i],
                                 genome = genome, length = 500,
                                 column_name_prefix = FALSE)

summit500 <- dplyr::select(summit500, peakId, summitSeq) %>% 
  dplyr::rename(summitSeq500 = summitSeq)


peakDf <- import_peaks_as_df(file = exptData$peakFile[i],
                             sampleId = exptData$sampleId[i],
                             peakFormat = "narrowPeak", rename = F)

peakDf <- dplyr::left_join(x = peakDf, y = summit200, by = "peakId") %>% 
  dplyr::left_join(y = summit500, by = "peakId")

readr::write_tsv(x = peakDf, path = paste(outPrefix, ".", exptData$sampleId[i], ".tab", sep = ""))















