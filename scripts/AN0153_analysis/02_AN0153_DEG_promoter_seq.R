suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(BSgenome.Anidulans.FGSCA4.AspGD))
suppressPackageStartupMessages(library(here))


## extract promoter sequence for AN0153_OE/WT polII DEGs

rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

##################################################################################

analysisName <- "AN0153_OE_vs_WT.DEG_promoter"
outDir <- here("analysis", "04_AN0153_analysis", "01_motif_enrichment")

outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
diffDataPath <- here::here("analysis", "06_polII_diff")

polIIDiffIds <- "AN0153_sCopy_OE_vs_MH11036"


if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

orgDb <- org.Anidulans.FGSCA4.eg.db
genome <- BSgenome.Anidulans.FGSCA4.AspGD
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- -1 * cutoff_lfc

col_lfc <- "log2FoldChange"
col_fdr <- "padj"

##################################################################################

## extract sequence for all gene promoter
geneInfo <- suppressMessages(
  AnnotationDbi::select(
    x = txDb, keys = AnnotationDbi::keys(x = txDb),
    columns = c("GENEID", "TXTYPE"), keytype = "GENEID")) %>%
  dplyr::rename(geneId = GENEID, txType = TXTYPE)

geneInfo <- dplyr::filter(geneInfo, !txType %in% c("tRNA", "rRNA", "snRNA", "snoRNA")) %>% 
  dplyr::filter(!grepl(pattern = "uORF", x = geneId))


genesGr <- GenomicFeatures::genes(x = txDb, filter = list(gene_id = geneInfo$geneId))

geneProGr <- GenomicRanges::promoters(x = genesGr, upstream = 500, downstream = 100)

geneProGr <- GenomicRanges::trim(geneProGr)

geneProGr <- geneProGr[which(width(geneProGr) == 600)]

geneProSeq <- Biostrings::getSeq(x = genome, names = geneProGr)

Biostrings::writeXStringSet(
  x = geneProSeq, format = "fasta",
  filepath = paste(outDir, "/A_nidulans_promoters.500_TSS_100.fasta", sep = "")
)


##################################################################################
## extract sequence for DEG promoters
rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison == polIIDiffIds)


diffData <- suppressMessages(
  readr::read_tsv(file = rnaseqInfo$deg)
) %>% 
  dplyr::filter(abs(!!sym(col_lfc)) > cutoff_lfc & !!sym(col_fdr) <= cutoff_fdr)


degGr <- GenomicFeatures::genes(x = txDb, filter = list(gene_id = diffData$geneId))

degPromoterSeq <- geneProSeq[intersect(diffData$geneId, names(geneProSeq))]

Biostrings::writeXStringSet(
  x = degPromoterSeq, format = "fasta",
  filepath = paste(outPrefix, ".fasta", sep = "")
)










