suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggridges))

## rank of SMTF in OE and WT polII-ChIPseq data: geom_point
## rank of SMTF in public RNAseq data from SRA: geom_density_ridges

rm(list = ls())
source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

##################################################################################
analysisName <- "SMTF_rank_SRA"
outDir <- here::here("analysis", "08_polII_analysis", "01_polII_DEGs_summary")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

polII_dataPath <- here::here("data", "polII_data")
diffDataPath <- here::here("analysis", "06_polII_diff")

file_genes <- here::here("data", "reference_data", "AN_genes_for_polII.bed")
file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_degIds <- here::here("data", "reference_data", "production_data.polII_DEG_ids.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
file_sraFpkm <- here::here("data", "Aspergillus_nidulans_FGSC_A4.RNAseq_SRA.FPKM.txt")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

##################################################################################

geneSet <- suppressMessages(readr::read_tsv(
  file = file_genes,
  col_names = c("chr", "start", "end", "geneId", "score", "strand"))) %>%
  dplyr::select(geneId)

degIds <- suppressMessages(readr::read_tsv(file = file_degIds))

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% degIds$degId)

polIIsamples <- unique(unlist(stringr::str_split(string = rnaseqInfo$samples, pattern = ";")))

polIIInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = polIIsamples,
  dataPath = polII_dataPath)

polIICols <- sapply(
  X = c("is_expressed", "rank", "perRank"),
  FUN = function(x){ structure(paste(x, ".", polIIInfo$sampleId, sep = ""), names = polIIInfo$sampleId) },
  simplify = F, USE.NAMES = T) %>% 
  append(
    list(exp = structure(polIIInfo$sampleId, names = polIIInfo$sampleId))
  )

##################################################################################
## prepare the data
polIIData <- get_polII_expressions(genesDf = geneSet, exptInfo = polIIInfo)

polIIRanks <- dplyr::select(polIIData, geneId, !!!polIIsamples) %>% 
  dplyr::mutate_at(
    .vars = vars(!matches("geneId")),
    .funs = ~rank(., ties.method = "min")
  )

smtf_polIIRanks <- dplyr::select(rnaseqInfo, geneId = SM_TF) %>% 
  dplyr::left_join(y = polIIRanks, by = "geneId") %>% 
  tidyr::pivot_longer(
    cols = -geneId, names_to = "sampleId", values_to = "rank"
  )

plotData <- dplyr::left_join(
  x = smtf_polIIRanks,
  y = dplyr::select(polIIInfo, sampleId, SM_TF, copyNumber, condition, timePoint),
  by = "sampleId"
) %>% 
  dplyr::filter(geneId == SM_TF | condition == "WT") %>% 
  tidyr::replace_na(replace = list(copyNumber = "WT")) %>% 
  dplyr::mutate(
    copyNumber = factor(copyNumber, levels = c("sCopy_OE", "mCopy_OE", "WT"))
  ) %>% 
  dplyr::arrange(geneId, copyNumber)


plotData$geneName <- AnnotationDbi::mapIds(
  x = orgDb, keys = plotData$geneId, column = "GENE_NAME", keytype = "GID"
)

plotData <- dplyr::mutate(
  plotData,
  geneLabel = paste(geneId, " (", geneName, ")", sep = ""),
  geneLabel = if_else(condition = geneId == geneName, true = geneId, false = geneLabel)
)

## add SRA FPKM data
sraData <- suppressMessages(readr::read_tsv(file = file_sraFpkm)) %>% 
  dplyr::rename(geneId = geneName)

sraRanks <- dplyr::mutate_at(
  .tbl = sraData, .vars = vars(!matches("geneId")),
  .funs = ~rank(., ties.method = "min")
) 

smtf_sraRanks <- dplyr::select(rnaseqInfo, geneId = SM_TF) %>% 
  dplyr::left_join(y = sraRanks, by = "geneId") %>% 
  tidyr::pivot_longer(
    cols = -geneId, names_to = "sampleId", values_to = "rank"
  )

smtf_sraRanks$geneName <- AnnotationDbi::mapIds(
  x = orgDb, keys = smtf_sraRanks$geneId, column = "GENE_NAME", keytype = "GID"
)

smtf_sraRanks <- dplyr::mutate(
  smtf_sraRanks,
  geneLabel = paste(geneId, " (", geneName, ")", sep = ""),
  geneLabel = if_else(condition = geneId == geneName, true = geneId, false = geneLabel)
)

sraRankMedian <- dplyr::group_by(smtf_sraRanks, geneName) %>% 
  dplyr::summarise(rankMedian = median(rank), rankMean = mean(rank)) %>% 
  dplyr::arrange(rankMedian)

plotData <- dplyr::mutate(
  plotData,
  geneName = forcats::fct_relevel(geneName, !!!sraRankMedian$geneName)
)

smtf_sraRanks <- dplyr::mutate(
  smtf_sraRanks,
  geneName = forcats::fct_relevel(geneName, !!!sraRankMedian$geneName)
)

##################################################################################
## plot the data
pt_oe_gene_rank <- ggplot(
  mapping = aes(x = rank, y = geneName)
) +
  geom_density_ridges(
    data = smtf_sraRanks
  ) +
  geom_point(data = plotData, mapping = aes(fill = copyNumber), size = 3, shape = 21) +
  geom_vline(xintercept = nrow(geneSet), linetype = "dashed") +
  coord_cartesian(xlim = c(0, nrow(geneSet))) +
  scale_x_continuous(
    breaks = c(1, 3000, 6000, 9000, nrow(geneSet)),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_fill_manual(
    values = c("sCopy_OE" = "red", "mCopy_OE" = "orange", "WT" = "green")
  ) +
  labs(
    title = "Rank of SM cluster TF gene in its own OE and WT polII ChIPseq data (points) and rank distribution of SMTF in public RNAseq data",
    # subtitle = paste("min rank = 0, max rank =", nrow(geneSet)),
    x = "rank(FPKM)"
  ) +
  guides(fill = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.title.y = element_blank()
  )


ggsave(filename = paste(outPrefix, ".png", sep = ""), plot = pt_oe_gene_rank, width = 14, height = 8)
ggsave(filename = paste(outPrefix, ".pdf", sep = ""), plot = pt_oe_gene_rank, width = 14, height = 8)








