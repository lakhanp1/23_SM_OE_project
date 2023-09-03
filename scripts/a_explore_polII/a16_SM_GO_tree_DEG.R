suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(GO.db))
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(Rgraphviz))

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/04_GO_enrichment/s01_enrichment_functions.R")
source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################
analysisName <- "temp_explore"
outDir <- here::here("analysis", "08_polII_analysis", "SM_GO_summary")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

if(!dir.exists(outDir)){
  dir.create(outDir, recursive = T)
}

file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
file_sampleInfo <- here::here("data", "reference_data", "polII_sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
diffDataPath <- here::here("analysis", "06_polII_diff")

cutoff_fpkm <- 10
cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- -1 * cutoff_lfc

col_lfc <- "log2FoldChange"

orgDb <- org.Anidulans.FGSCA4.eg.db
keggOrg <- 'ani'
col_degId <- "geneId"
col_degOrgdbKey <- "GID"
col_kegg <- "KEGG_ID"
col_gsea <- "geneId"
col_geneName <- "GENE_NAME"

###########################################################################

smGoId <- c("GO:0019748")

levelNodes <- AnnotationDbi::mget(x = c(smGoId), envir = GO.db::GOBPOFFSPRING, ifnotfound = NA) %>% 
  unlist() %>% 
  unique()

## select only GO terms which are annotated in org.db
levelNodes <- intersect(levelNodes, keys(x = orgDb, keytype = "GOALL"))

levelNodes <- setdiff(x = levelNodes, y = c("GO:0019748"))

geneToGo <- AnnotationDbi::select(
  x = orgDb, keys = levelNodes, columns = "GID", keytype = "GOALL"
)

length(unique(geneToGo$GID))

tmpGoList <- AnnotationDbi::mapIds(
  x = orgDb, keys = c(smGoId, "GO:0044550"), keytype = "GOALL", column = "GID", multiVals = list
)

goGeneCount <- dplyr::group_by(geneToGo, GOALL) %>% 
  dplyr::summarise(nGenes = n_distinct(GID))

## build a GO table 
goTable <- suppressMessages(
  AnnotationDbi::select(
    x = GO.db,
    keys = levelNodes,
    columns = c("GOID", "ONTOLOGY", "TERM"),
    keytype = "GOID")) %>% 
  dplyr::left_join(y = goGeneCount, by = c("GOID" ="GOALL"))

dplyr::arrange(goTable, desc(nGenes)) %>% 
  readr::write_tsv(file = paste(outPrefix, ".SM_GO_stats.tab", sep = ""))


parentNodes <- AnnotationDbi::mget(
  x = goTable$GOID, envir = GO.db::GOBPPARENTS
)

## use org.db object to extract gene->GO list
geneID2GO <- suppressMessages(
  AnnotationDbi::mapIds(
    x = orgDb, keys = keys(x = orgDb, keytype = col_degOrgdbKey),
    column = "GOALL", keytype = col_degOrgdbKey,
    multiVals = list
  )
)

geneNames <- names(geneID2GO)
genes <- unname(
  unlist(
    AnnotationDbi::mapIds(
      x = orgDb, keys = smGoId, column = "GID", keytype = "GOALL", multiVals = list
    )
  ))

geneList <- factor(as.integer(geneNames %in% genes))
names(geneList) <- geneNames

goData <- suppressMessages(
  new(Class = "topGOdata",
      ontology = "BP",
      allGenes = geneList,
      annot = annFUN.gene2GO,
      gene2GO = geneID2GO,
      nodeSize = 1)
)

# resFisherWeight <- suppressMessages(
#   topGO::runTest(goData, algorithm = "weight01", statistic = "fisher")
# )
# 
# 
# 
# ## get the result table, filter by pValue cutoff 0.05 and calculate Rich factor
# resultTab <- topGO_enrichment(
#   orgdb = orgDb, genes = genes, type = "BP", inKeytype = "GID"
# )
# 
# topGO::printGraph(
#   object = goData, result = resFisherWeight, firstSigNodes = 5, 
#   fn.prefix = outPrefix, pdfSW = TRUE
# )
# 
# showSigOfNodes(
#   goData, score(resFisherWeight), useInfo = "all"
# )


dag <- subGraph(snodes = c(levelNodes, smGoId), graph = topGO::graph(goData))
plot(dag)

dag2 <- graph::addNode(node = unique(geneToGo$GID), object = dag)
dag2 <- graph::addEdge(from = geneToGo$GID, to = geneToGo$GOALL, graph = dag2)

plot(dag2)

graph::degree(object = dag2, Nodes = c(levelNodes, smGoId))

smGenes <- AnnotationDbi::select(
  x = orgDb, keys = keys(x = orgDb, keytype = "SM_ID"), columns = "GID", keytype = "SM_ID"
)














