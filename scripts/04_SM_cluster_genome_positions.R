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

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
file_smBed <- "E:/Chris_UM/Database/A_Nidulans/SM_genes.bed"

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

##################################################################################
## process the data
productionData <- suppressMessages(readr::read_tsv(file = file_productionData))

prodDataFiltered <- productionData %>% 
  dplyr::filter(has_polII_ChIP == "has_data", has_TF_ChIP == "has_data", copyNumber == "sCopy")


smInfo <- suppressMessages(
  AnnotationDbi::select(
    x = orgDb, keys = keys(orgDb, keytype = "SM_ID"),
    columns = c("GID", "SM_ID", "SM_CLUSTER"), keytype = "SM_ID")
) %>% 
  dplyr::rename(geneId = GID) %>% 
  dplyr::mutate(hasData = "no_TF")


## BGCs with TF
bgcWithTf <- productionData %>% 
  dplyr::filter(!is.na(SM_ID)) %>%
  dplyr::select(SM_ID) %>% 
  dplyr::distinct() %>% 
  tibble::deframe()

## BGCs with TF+polII data
bgcWithData <- productionData %>% 
  dplyr::filter(
    has_polII_ChIP == "has_data", has_TF_ChIP == "has_data",
    copyNumber == "sCopy", !is.na(SM_ID)
  ) %>%
  dplyr::select(SM_ID) %>% 
  dplyr::distinct() %>% 
  tibble::deframe()

smInfo$hasData[which(smInfo$SM_ID %in% bgcWithTf)] <- "has_TF"
smInfo$hasData[which(smInfo$SM_ID %in% bgcWithData)] <- "has_TF_data"

smGr <- GenomicFeatures::genes(
  x = txDb, filter = list(gene_id = smInfo$geneId)
) %>% 
  as.data.frame %>% 
  dplyr::left_join(y = smInfo, by = c("gene_id" = "geneId")) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

seqlevels(smGr) <- seqlevels(txDb)
seqinfo(smGr) <- seqinfo(txDb)


## extract individual cluster range as GRanges and save as BED
clusterGrl <- GenomicRanges::split(x = smGr, f = smGr$SM_ID)
clusterGr <- sort(unlist(range(clusterGrl, ignore.strand = TRUE)))
mcols(clusterGr)$SM_ID <- names(clusterGr)

mcols(clusterGr)$SM_CLUSTER <- AnnotationDbi::mapIds(
  x = orgDb, keys = mcols(clusterGr)$SM_ID,
  column = "SM_CLUSTER", keytype = "SM_ID"
)

mcols(clusterGr)$name <- paste(
  mcols(clusterGr)$SM_ID, "::", mcols(clusterGr)$SM_CLUSTER, sep = ""
) %>% 
  stringr::str_replace_all(pattern = "\\s+", replacement = "_")

rtracklayer::export.bed(
  object = clusterGr,
  con = here::here("data", "reference_data", "SM_cluster_regions.bed")
)

##################################################################################
## combine SM genes into regions
regionsGr <- sort(reduce(clusterGr), ignore.strand = TRUE)
mcols(regionsGr)$regionId <- paste(
  "region_", sprintf(fmt = "%02d",1:length(regionsGr)), sep = ""
)

mcols(regionsGr)$region <- as.character(regionsGr, ignore.strand = TRUE)

genePos <- sort(
  GenomicFeatures::genes(x = txDb, filter = list(gene_id = unique(smInfo$geneId))),
  ignore.strand = TRUE
)

ovlp_genes <- GenomicRanges::findOverlaps(
  query = genePos, subject = regionsGr, type = "within"
)

ovlp_clusters <- GenomicRanges::findOverlaps(
  query = clusterGr, subject = regionsGr
)

genesToCluster <- dplyr::select(smInfo, geneId, SM_ID) %>% 
  dplyr::group_by(geneId) %>% 
  dplyr::summarise(
    SM_IDs = paste(SM_ID, collapse = ";"),
    nClusters = n()
  )

mcols(genePos)$regionId[ovlp_genes@from] <- mcols(regionsGr)$regionId[ovlp_genes@to]
mcols(genePos)$region <- mcols(regionsGr)$region[ovlp_genes@to]

geneRegions <- as.data.frame(mcols(genePos)) %>% 
  dplyr::rename(geneId = gene_id) %>% 
  dplyr::left_join(y = genesToCluster, by = "geneId")

readr::write_tsv(
  x = geneRegions,
  file = here::here("data", "reference_data", "SM_genes_regions.tab")
)

##################################################################################
## plot the data
pt_theme <- theme(
  plot.title = element_text(size = 18, face = "bold", hjust = 1),
  axis.text = element_text(size = 16, face = "bold", color = "black"),
  axis.ticks.length = unit(2, "mm"),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title = element_blank(),
  legend.position = "bottom",
  legend.text = element_text(size = 16, face = "bold"),
  strip.background = element_rect(fill = "white"),
  strip.text.y = element_text(angle = 0),
  panel.border = element_blank(),
  panel.grid = element_blank()
)

pt_all <- autoplot(object = seqinfo(smGr), layout = "karyogram") +
  ggbio::layout_karyogram(
    data = smGr, geom = "rect", color = "black", 
    mapping = aes(fill = "black"),
    ylim = c(2, 8)
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("black" = "black"),
    labels = c("black" = "BGCs")
  ) +
  labs(title = "Distribution of BGCs in A. nidulans genome") +
  theme_clear() +
  pt_theme

pdf(file = paste(outPrefix, ".pdf", sep = ""), width = 12, height = 6)
pt_all
dev.off()



pt_data <- autoplot(object = seqinfo(smGr), layout = "karyogram") +
  ggbio::layout_karyogram(
    data = smGr, geom = "rect", 
    mapping = aes(color = hasData, fill = hasData),
    ylim = c(2, 8)
  ) +
  scale_color_manual(
    name = NULL,
    values = c("has_TF_data" = "blue", "has_TF" = "black", "no_TF" = "grey"),
    labels = c("has_TF_data" = "TF + data", "has_TF" = "TF & no data", "no_TF" = "no TF")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("has_TF_data" = "blue", "has_TF" = "black", "no_TF" = "grey"),
    labels = c("has_TF_data" = "TF + data", "has_TF" = "TF & no data", "no_TF" = "no TF")
  ) +
  labs(title = "Distribution of BGCs in A. nidulans genome") +
  theme_clear() +
  pt_theme


pdf(file = paste(outPrefix, ".data.pdf", sep = ""), width = 12, height = 6)
pt_data
dev.off()




