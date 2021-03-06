suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggrepel))


rm(list = ls())


##################################################################################
analysisName <- "OE_TF_polII_rank"
outDir <- here::here("analysis", "08_polII_analysis", "01_polII_DEGs_summary")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")

file_genes <- here::here("data", "reference_data", "AN_genes_for_polII.bed")
orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

polII_dataPath <- here::here("data", "polII_data")

file_polII <- paste(polII_dataPath, "/", "samples_polII.DESeq2.list", sep = "")

##################################################################################

geneSet <- suppressMessages(readr::read_tsv(
  file = file_genes,
  col_names = c("chr", "start", "end", "geneId", "score", "strand"))) %>%
  dplyr::select(geneId)

polIIsamples <- readr::read_tsv(file = file_polII, comment = "#")

polIIInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = polIIsamples$sampleId,
  dataPath = polII_dataPath)

polIICols <- sapply(
  X = c("is_expressed", "rank", "perRank"),
  FUN = function(x){ structure(paste(x, ".", polIIInfo$sampleId, sep = ""), names = polIIInfo$sampleId) },
  simplify = F, USE.NAMES = T) %>% 
  append(
    list(exp = structure(polIIInfo$sampleId, names = polIIInfo$sampleId))
  )

##################################################################################


oeGeneRanks <- NULL

i <- 1

## get polII signal rank of TF gene in its own OE dataset
for (i in 1:nrow(polIIInfo)) {
  
  if(is.na(polIIInfo$SM_TF[i])){
    next
  }
  
  dt <- get_polII_expressions(genesDf = geneSet, exptInfo = polIIInfo[i, ]) %>% 
    dplyr::mutate(
      rank = rank(!!sym(polIIInfo$sampleId[i])),
      perRank = percent_rank(!!sym(polIIInfo$sampleId[i]))
    ) %>% 
    dplyr::filter(geneId == !!polIIInfo$SM_TF[i])
  
  
  oeGeneRanks <- dplyr::bind_rows(
    oeGeneRanks,
    dplyr::select(polIIInfo[i, ], sampleId, SM_TF, copyNumber, condition, timePoint) %>% 
      dplyr::mutate(
        geneId = dt$geneId,
        polII_signal = dt[[polIIInfo$sampleId[i]]],
        rank = dt$rank,
        perRank = dt$perRank)
  )
  
  
}


## add TF gene rank in WT polII data
wtPolIIInfo <- dplyr::filter(polIIInfo, condition == "WT" & timePoint == "16h")

wtPolIIData <- get_polII_expressions(genesDf = geneSet, exptInfo = wtPolIIInfo)

rankFormula <- purrr::map(
  .x = wtPolIIInfo$sampleId,
  .f = function(x){
    
    nameRank <- sym(paste("rank.", x, sep = ""))
    namePerRank <- sym(paste("perRank.", x, sep = ""))
    
    list(
      quos(!!nameRank := rank(!!sym(x))),
      quos(!!namePerRank := percent_rank(!!sym(x)))
    )
  }
) %>% 
  unlist()

wtPolIIData <- dplyr::mutate(wtPolIIData, !!!rankFormula)

wtRanks <- data.table::as.data.table(wtPolIIData) %>% 
  data.table::melt(
    measure.vars = purrr::map(.x = polIICols, .f = ~ unname(.x[wtPolIIInfo$sampleId])),
    variable.name = "sampleId"
  ) %>% 
  as_tibble() %>% 
  dplyr::select(-is_expressed, polII_signal = exp)

wtRanks$sampleId <- forcats::fct_recode(
  wtRanks$sampleId,
  !!!structure(levels(wtRanks$sampleId), names = wtPolIIInfo$sampleId)
)


wtGeneRanks <- dplyr::left_join(x = wtRanks, y = wtPolIIInfo, by = "sampleId") %>% 
  dplyr::select(!!!colnames(oeGeneRanks)) %>% 
  tidyr::replace_na(replace = list(copyNumber = "WT")) %>% 
  dplyr::filter(geneId %in% oeGeneRanks$geneId)


plotData <- dplyr::bind_rows(oeGeneRanks, wtGeneRanks) %>% 
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

## plot the data
pt_oe_gene_rank <- dplyr::filter(plotData, copyNumber != "deletion") %>% 
  ggplot(mapping = aes(x = geneLabel, y = rank, color = copyNumber)) +
  geom_hline(yintercept = nrow(geneSet), linetype = "dashed") +
  geom_point(size = 3) +
  coord_cartesian(ylim = c(0, nrow(geneSet))) +
  scale_y_continuous(
    breaks = c(1, 3000, 6000, 9000, nrow(geneSet)),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_color_manual(
    values = c("sCopy_OE" = "red", "mCopy_OE" = "#4daf4a", "WT" = "black")
  ) +
  labs(
    title = "Rank of SM cluster TF gene in its own over-expression and WT polII ChIPseq data",
    subtitle = paste("min rank = 0, max rank =", nrow(geneSet)),
    x = "SM TF", y = "rank(polII ChIPseq signal)"
  ) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title = element_text(size = 14)
  )


ggsave(filename = paste(outPrefix, ".png", sep = ""), plot = pt_oe_gene_rank, width = 12, height = 8)
ggsave(filename = paste(outPrefix, ".pdf", sep = ""), plot = pt_oe_gene_rank, width = 12, height = 8)














