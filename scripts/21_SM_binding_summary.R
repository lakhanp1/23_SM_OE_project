library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(TxDb.Anidulans.FGSCA4.AspGD.GFF)
library(here)

## TF binding summary plot at SM gene clusters

rm(list = ls())


##################################################################################
analysisName <- "SM_cluster_summary"
outDir <- here::here("analysis", "06_SM_cluster_summary")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_exptInfo <- here::here("data", "referenceData/sample_info.txt")

file_genes <- here::here("data", "referenceData/AN_genes_for_polII.bed")
orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

TF_dataPath <- here::here("data", "TF_data")

file_tf_macs2 <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")
file_tf <- paste(TF_dataPath, "/", "sample_tf.list", sep = "")


smInfo <- suppressMessages(
  AnnotationDbi::select(x = orgDb, keys = keys(orgDb, keytype = "SM_CLUSTER"),
                        columns = c("GID", "SM_ID"), keytype = "SM_CLUSTER")) %>% 
  dplyr::rename(geneId = GID)


genePosition <- GenomicFeatures::genes(x = txDb, filter = list(gene_id = unique(smInfo$geneId))) %>% 
  as.data.frame() %>% 
  dplyr::rename(chr = seqnames) %>% 
  dplyr::select(-width)

smInfo <- dplyr::left_join(x = smInfo, y = genePosition, by = c("geneId" = "gene_id")) %>% 
  dplyr::arrange(chr, start) %>% 
  dplyr::mutate(
    rowId = row_number(),
    geneId = forcats::as_factor(geneId),
    SM_CLUSTER = forcats::as_factor(SM_CLUSTER)
  )

# smInfo <- smInfo[1:56, ]
# dplyr::group_by(geneId) %>%
# tidyr::nest()
# dplyr::mutate(SM_CLUSTER = paste(SM_CLUSTER, collapse = ";"),
#               SM_ID = paste(SM_ID, collapse = ";")) %>% 
# dplyr::slice(1L) %>% 
# dplyr::ungroup()

##################################################################################


tfSampleList <- readr::read_tsv(file = file_tf_macs2, col_names = c("sampleId"),  comment = "#")

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSampleList$sampleId,
  dataPath = TF_dataPath)

tfInfoList <- purrr::transpose(tfInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))


tfCols <- sapply(
  c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
    "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
    "peakCoverage", "pvalFiltered", "summitSeq"),
  FUN = function(x){ structure(paste(x, ".", tfInfo$sampleId, sep = ""), names = tfInfo$sampleId) },
  simplify = F, USE.NAMES = T)


i <- 1
plotData <- NULL

for (i in 1:nrow(tfInfo)) {
  
  cat(i, ":", tfInfo$sampleId[i], "\n")
  
  if(file.exists(tfInfo$peakAnno[i])){
    
    peakTargets <- suppressMessages(
      readr::read_tsv(file = tfInfo$peakAnno[i])
    ) %>% 
      dplyr::filter(peakPval >= 20) %>% 
      dplyr::mutate(hasPeak = TRUE) %>% 
      dplyr::group_by(geneId) %>% 
      dplyr::arrange(desc(peakPval), .by_group = TRUE) %>% 
      dplyr::slice(1L) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(geneId, peakPval, peakType, peakPosition, hasPeak)
    
    
    smClusterTargets <- dplyr::left_join(x = smInfo, y = peakTargets, by = "geneId") %>% 
      dplyr::mutate(
        sampleId = tfInfo$sampleId[i],
        sm_tf = tfInfo$SM_TF[i],
        timepoint = tfInfo$timePoint[i],
        condition = tfInfo$condition[i]
      ) %>% 
      tidyr::replace_na(replace = list(hasPeak = FALSE))
    
    plotData <- dplyr::bind_rows(plotData, smClusterTargets)
  }
  
}



plotData <- plotData %>% dplyr::mutate(
  geneId = forcats::as_factor(geneId),
  sampleId = forcats::as_factor(sampleId)
)


pt_binding <- ggplot(data = plotData, mapping = aes(x = geneId, y = sampleId)) +
  # geom_point(mapping = aes(fill = hasPeak), shape = 24, color = "white", size = 2) +
  geom_tile(mapping = aes(fill = hasPeak), height = 0.8) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = alpha("white", 0)), guide = FALSE) +
  facet_grid(cols = vars(SM_CLUSTER), scales = "free", space = "free_x") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_blank()
  )


png(filename = paste(outPrefix, ".binding.png", sep = ""), width = 10000, height = 4000, res = 120)

pt_binding

dev.off()











