library(chipmine)
library(org.Anidulans.eg.db)
library(TxDb.Anidulans.AspGD.GFF)
library(here)
library(ggridges)


## 1) peak enrichment distribution

rm(list = ls())

##################################################################################

file_exptInfo <- here::here("data", "referenceData/sample_info.txt")

file_genes <- here::here("data", "referenceData/AN_genes_for_polII.bed")
orgDb <- org.Anidulans.eg.db
txDb <- TxDb.Anidulans.AspGD.GFF

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")
other_dataPath <- here::here("data", "other_data")

file_tf_macs2 <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")
file_tf <- paste(TF_dataPath, "/", "sample_tf.list", sep = "")

matrixType <- "normalizedmatrix"
matrixDim = c(200, 200, 100, 10)

geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "geneId", "score", "strand")) %>%
  dplyr::select(-score) %>%
  dplyr::mutate(length = end - start)

geneSet <- GenomicFeatures::genes(x = txDb, columns = c("gene_id", "tx_id", "tx_name"),
                                  filter = list(gene_id = geneSet$geneId)) %>% 
  as.data.frame() %>% 
  dplyr::select(geneId = gene_id, chr = seqnames, start, end, strand)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$geneId, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID"))


smGenes <- suppressMessages(
  AnnotationDbi::select(
    x = orgDb,
    keys = keys(x = orgDb, keytype = "SM_GENE"),
    columns = c("SM_GENE", "SM_CLUSTER", "SM_ID"),
    keytype = "SM_GENE")) %>% 
  dplyr::mutate(geneType = "SM")


##################################################################################

tfSampleList <- readr::read_tsv(file = file_tf_macs2, col_names = c("id"),  comment = "#")

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSampleList$id,
  dataPath = TF_dataPath,
  profileMatrixSuffix = matrixType)

tfInfoList <- purrr::transpose(tfInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))


i <- 2

backboneGene <- tfInfo$SM_TF[i]

## few TFs are mapped to multiple SM clusters. So preparing the list of SM tf data
smTfInfo <- suppressMessages(
  AnnotationDbi::select(
    x = orgDb,
    keys = backboneGene,
    columns = c("SM_GENE", "SM_CLUSTER", "SM_ID"),
    keytype = "SM_GENE")) %>% 
  dplyr::filter(!is.na(SM_ID))


smGeneCategory <- smGenes %>% 
  dplyr::mutate(geneType = if_else(condition = SM_ID %in% !!smTfInfo$SM_ID,
                                   true = "cluster-SM", false = geneType))


peakAnno <- import_peak_annotation(sampleId = tfInfo$sampleId[i],
                                   peakAnnoFile = tfInfo$peakAnno[i],
                                   renameColumn = FALSE)

# intersect(peakAnno$geneId, smGenes$SM_GENE)

peakAnno <- dplyr::left_join(x = peakAnno, y = smGeneCategory, by = c("geneId" = "SM_GENE")) %>% 
  tidyr::replace_na(replace = list(geneType = "non-SM")) %>% 
  dplyr::mutate(
    drawOrder = dplyr::case_when(
      geneType == "cluster-SM" ~ 1,
      geneType == "SM" ~ 2,
      geneType == "non-SM" ~ 3
    )
  ) %>% 
  dplyr::distinct()

peakAnno$geneType <- factor(x = peakAnno$geneType, levels = c("non-SM", "SM", "cluster-SM"))

summary(peakAnno$peakPval)
summary(peakAnno$peakEnrichment)

# hist(peakAnno$peakEnrichment,
#      breaks = quantile(peakAnno$peakEnrichment, seq(from = 0, to = 1, by = 0.1)))


ggplot(data = peakAnno,
       mapping = aes(x = tfInfo$sampleId[i], y = peakEnrichment, fill = geneType)) +
  geom_violin(fill = "green") +
  geom_jitter(mapping = aes(fill = geneType, alpha = geneType),
              shape=21, size = 2,
              position=position_jitter(width = 0.4, seed = 1)) +
  # geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.6,
  #              binpositions = "all", binwidth = 0.2, stackratio=1) +
  scale_fill_manual(name = "peak type",
                    values = c("cluster-SM" = "red", "SM" = "blue", "non-SM" = "grey")) +
  scale_alpha_manual(values = c("cluster-SM" = 1, "SM" = 0.6, "non-SM" = 0.4)) +
  guides(alpha = FALSE) +
  theme_bw()


ggplot(data = peakAnno,
       mapping = aes(x = peakEnrichment, y = tfInfo$sampleId[i])) +
  geom_density_ridges(
    # mapping = aes(point_fill = geneType, fill = geneType),
    alpha = 0.5,
    jittered_points = TRUE,
    # position = "raincloud",
    point_shape = 21,
    point_color = "black", point_alpha = 0.5
  )









