library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(TxDb.Anidulans.FGSCA4.AspGD.GFF)
library(BSgenome.Anidulans.FGSCA4.AspGD)
library(here)


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

matrixType <- "2kb_summit"
up <- 2000
down <- 2000
body <- 0
bin <- 10
matrixDim = c(c(up, body, down)/bin, bin)

##################################################################################

geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "geneId", "score", "strand")) %>%
  dplyr::select(geneId)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$geneId, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID"))


txInfo <- suppressMessages(
  AnnotationDbi::select(
    x = txDb, keys = AnnotationDbi::keys(x = txDb, keytype = "TXID"),
    columns = c("GENEID", "TXNAME", "TXTYPE"), keytype = "TXID")) %>%
  dplyr::mutate(TXID = as.character(TXID)) %>%
  dplyr::rename(geneId = GENEID, txName = TXNAME, txType = TXTYPE)

txInfo <- dplyr::filter(txInfo, !txType %in% c("tRNA", "rRNA", "snRNA", "snoRNA")) %>% 
  dplyr::filter(!grepl(pattern = "uORF", x = geneId))


tfSampleList <- readr::read_tsv(file = file_tf_macs2, col_names = c("sampleId"),  comment = "#")

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSampleList$sampleId,
  dataPath = TF_dataPath,
  profileMatrixSuffix = matrixType)

##################################################################################

i <- 3

for(i in 1:nrow(tfInfo)){
  
  cat(i, ":", tfInfo$sampleId[i], "\n")
  
  ## annotate peaks and prepare gene level annotation file
  peakType <- dplyr::case_when(
    tfInfo$peakType[i] == "narrow" ~ "narrowPeak",
    tfInfo$peakType[i] == "broad" ~ "broadPeak"
  )
  
  peakAn <- narrowPeak_annotate(
    peakFile = tfInfo$peakFile[i],
    txdb = txDb,
    txIds = txInfo$TXID,
    fileFormat = peakType,
    promoterLength = 1000,
    upstreamLimit = 1500,
    bidirectionalDistance = 1000,
    includeFractionCut = 0.7,
    bindingInGene = FALSE,
    insideSkewToEndCut = 0.7,
    removePseudo = TRUE,
    output = tfInfo$peakAnno[i])
  
  if( !is.null(peakAn) ){
    tfDf <- gene_level_peak_annotation(
      sampleId = tfInfo$sampleId[i],
      peakAnnotation = tfInfo$peakAnno[i],
      genesDf = geneSet,
      peakFile = tfInfo$peakFile[i],
      bwFile = tfInfo$bwFile[i],
      outFile = tfInfo$peakTargetFile[i])
  }
  
  
  # ## create profile matrix of 2kb region around peak summit
  # peaksGr <- rtracklayer::import(con = tfInfo$peakFile[i], format = peakType)
  # 
  # if(length(peaksGr) > 0){
  #   if(peakType == "broadPeak"){
  #     mcols(peaksGr)$peak <- round(width(peaksGr) / 2)
  #   }
  # 
  #   peakSummitGr <- GenomicRanges::narrow(x = peaksGr,
  #                                         start = pmax(peaksGr$peak, 1),
  #                                         width = 1)
  # 
  #   profileMat <- bigwig_profile_matrix(bwFile = tfInfo$bwFile[i],
  #                                       regions = peakSummitGr,
  #                                       signalName = tfInfo$profileName[i],
  #                                       extend = c(up, down),
  #                                       targetName = "summit",
  #                                       storeLocal = T,
  #                                       localPath = tfInfo$matFile[i])
  # }
  
}



