suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(viridis))

## intergenic region types and width distribution 
## SMTF peak overlap with intergenic regions

rm(list = ls())

##################################################################################

analysisName <- "intergenic_regions"
outDir <- here::here("analysis", "02_QC_TF", "intergenic_regions")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
file_exptInfo <- here::here("data", "reference_data/sample_info.txt")
file_genes <- here::here("data", "reference_data/AN_genes_for_polII.bed")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

TF_dataPath <- here::here("data", "TF_data")

cutoff_macs2Pval <- 20

matrixType <- "2kb_summit"
up <- 2000
down <- 2000
body <- 0
bin <- 10
matrixDim = c(c(up, body, down)/bin, bin)

wrap60 <- wrap_format(80)
##################################################################################

if(!dir.exists(outDir)){
  dir.create(path = outDir, recursive = TRUE)
}

productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
  dplyr::filter(has_polII_ChIP == "has_data", has_TF_ChIP == "has_data", copyNumber == "sCopy") %>% 
  dplyr::arrange(SM_ID, SMTF)

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = productionData$tfId,
  dataPath = TF_dataPath)

tfInfoList <- purrr::transpose(tfInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

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

## extract bidirectional regions
chrGr <- GenomicRanges::GRanges(
  seqnames = names(seqlengths(txDb)),
  ranges = IRanges(start = 1, width = seqlengths(txDb))
)

intergenicGr <- GenomicRanges::setdiff(
  x = chrGr, y = sort(unstrand(transcriptsGr))
)

mcols(intergenicGr)$name <- paste("intergenic_", 1:length(intergenicGr), sep = "")
mcols(intergenicGr)$leftStrand <- NA
mcols(intergenicGr)$rightStrand <- NA
mcols(intergenicGr)$leftGene <- "NA"
mcols(intergenicGr)$rightGene <- "NA"

gapLeftGr <- GenomicRanges::follow(
  x = intergenicGr,
  subject = transcriptsGr,
  ignore.strand = TRUE,
  select = "all"
)

gapRightGr <- GenomicRanges::precede(
  x = intergenicGr,
  subject = transcriptsGr,
  ignore.strand = TRUE,
  select = "all"
)


mcols(intergenicGr)$leftStrand[gapLeftGr@from] <- strand(transcriptsGr)[gapLeftGr@to]
mcols(intergenicGr)$rightStrand[gapRightGr@from] <- strand(transcriptsGr)[gapRightGr@to]

mcols(intergenicGr)$leftGene[gapLeftGr@from] <- transcriptsGr$gene_id[gapLeftGr@to]
mcols(intergenicGr)$rightGene[gapRightGr@from] <- transcriptsGr$gene_id[gapRightGr@to]

gapDf <- as.data.frame(intergenicGr) %>% 
  dplyr::select(-strand) %>% 
  dplyr::filter(!is.na(leftStrand), !is.na(rightStrand)) %>% 
  dplyr::mutate(
    gapType = dplyr::case_when(
      leftStrand == rightStrand ~ "tandem",
      # leftStrand == "+" & rightStrand == "+" ~ ">=(+)=>     >=(+)=>",
      # leftStrand == "-" & rightStrand == "-" ~ "<=(-)=<     <=(-)=<",
      leftStrand == "+" & rightStrand == "-" ~ "convergent",
      leftStrand == "-" & rightStrand == "+" ~ "divergent"
      # leftStrand == "+" & rightStrand == "-" ~ ">=(+)=>     <=(-)=<",
      # leftStrand == "-" & rightStrand == "+" ~ "<=(-)=<     >=(+)=>",
    )
  )

intergenicGr <- GenomicRanges::makeGRangesFromDataFrame(
  df = gapDf, keep.extra.columns = TRUE
)


itgSummary <- gapDf %>% 
  dplyr::group_by(gapType) %>% 
  dplyr::tally(name = "regionsN")


fun_summary <- function(x){
  return(
    data.frame(ymin = quantile(x, 0.25), y = quantile(x, 0.5), ymax = quantile(x, 0.75))
  )
}


## plot intergenic region width statistics
pt_width <- ggplot(
  data = gapDf,
  mapping = aes(x = gapType, y = pmin(width, 3000))
) +
  geom_quasirandom() +
  stat_summary(fun.data = fun_summary, color="red", size = 1) +
  scale_x_discrete(
    labels = c(
      "tandem" = "tandem",
      "convergent" = "convergent",
      "divergent" = "divergent"
    )
  ) +
  geom_text(
    data = itgSummary,
    mapping = aes(x = gapType, y = Inf, label = paste("n =", regionsN)),
    size = 6, vjust = 1
  ) +
  labs(
    x = "Gene orientation around intergenic region",
    y = "Intergenic region width",
    title = "Intergenic region width and orientation summary"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(2, "mm")
  )

ggplot2::ggsave(filename = paste(outPrefix, ".summary.pdf", sep = ""),
                plot = pt_width, width = 12, height = 8)

##################################################################################
## for each TF ChIPseq peaks, find the % of peaks in different intergenic regions
rowId <- 1

itgPeakStats <- NULL

for (rowId in 1:nrow(productionData)) {
  
  tfId <- productionData$tfId[rowId]
  peakType <- tfInfoList[[tfId]]$peakType
  
  peaksGr <- rtracklayer::import(con = tfInfoList[[tfId]]$peakFile, format = peakType)
  
  peaksGr <- peaksGr[peaksGr$pValue >= cutoff_macs2Pval, ]
  
  if(length(peaksGr) == 0){
    # cat(tfInfo$sampleId[rowId], "\n")
    next()
  }
  
  ## use peak summit instead of whole peak
  peakSummitGr <- GenomicRanges::resize(
    x = GenomicRanges::shift(x = peaksGr, shift = peaksGr$peak),
    width = 1, fix = "start"
  )
  
  
  ## peaks falling within intergenic regions
  ovHits <- GenomicRanges::findOverlaps(
    query = peakSummitGr,
    subject = intergenicGr,
    select = "all", ignore.strand = TRUE,
    type = "within"
  )
  
  
  # peakSummitGr[ovHits@from]
  ovlpGr <- as.data.frame(mcols(intergenicGr)[ovHits@to, ]) %>% 
    dplyr::mutate(
      peakId = peakSummitGr$name[ovHits@from]
    )
  
  ## summary stats of peaks in intergenic regions
  tfSummary <- ovlpGr %>% 
    dplyr::group_by(gapType) %>% 
    dplyr::tally(name = "intergenicPeaks")
  
  ## add a row for non-intergenic row
  tfSummary <- tfSummary %>% 
    dplyr::bind_rows(
      tibble(
        gapType = "non-intergenic",
        intergenicPeaks = length(peakSummitGr) - sum(tfSummary$intergenicPeaks)
      )
    ) %>% 
    dplyr::mutate(
      totalIntergenicPeaks = length(ovHits),
      totalPeaks = length(peakSummitGr),
      sampleId = tfInfo$sampleId[rowId],
      OESMTF = !!productionData$SMTF[rowId],
      OESMTF_name = !!productionData$SMTF_name[rowId]
    )
  
  itgPeakStats <- dplyr::bind_rows(itgPeakStats, tfSummary) 
}


itgPeakStats <- dplyr::left_join(
  x = itgPeakStats, y = itgSummary, by = "gapType"
)


itgPeakStats <- dplyr::mutate(
  itgPeakStats,
  fraction = round(x = intergenicPeaks/totalPeaks, digits = 3),
  gapType = forcats::fct_relevel(gapType, "tandem", "convergent", "divergent")
)

readr::write_tsv(x = itgPeakStats, file = paste(outPrefix, ".all_TF_distriution.data.tab", sep = ""))

pt2 <- ggplot(
  data = itgPeakStats,
  mapping = aes(x = gapType, y = fraction)
) +
  geom_quasirandom(mapping = aes(fill = totalPeaks), shape = 21, size = 3) +
  geom_text_repel(
    data = dplyr::filter(itgPeakStats, fraction >= 0.5),
    mapping = aes(x = gapType, y = fraction, label = OESMTF_name),
    size = 5, fontface = "italic",
    seed = 42, box.padding = 0.5,
    min.segment.length = 0,
    segment.size = unit(1, "mm")
  ) +
  stat_summary(fun.data = fun_summary, color="red", size = 1) +
  scale_fill_viridis(name = "log2(#peaks)", option = "viridis", direction = -1, trans = "log2") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Gene orientation around intergenic region",
    y = "% of peaks for a TF ChIPseq data",
    title = str_wrap("Distribution of TF ChIPseq datasets showing % of peaks in different chromosome regions"),
    subtitle = str_wrap(
      string = "Each point is one TF ChIPseq dataset. For each TF ChIPseq data, fraction of total peaks in each category of region is determined. These values are plotted as distribution for all the TF ChIPseq samples",
      width = 100)
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(2, "mm"),
    legend.text = element_text(size = 16),
    legend.key.height = unit(1, "cm"),
    legend.title = element_text(size = 16, face = "bold", angle = 90, vjust = 1)
  )

ggplot2::ggsave(filename = paste(outPrefix, ".all_TF_distriution.pdf", sep = ""),
                plot = pt2, width = 12, height = 8)



