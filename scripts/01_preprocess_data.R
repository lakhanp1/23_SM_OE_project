library(chipmine)
library(org.Anidulans.eg.db)
library(foreach)
library(doParallel)
library(here)

rm(list = ls())


# cl <- makeCluster(4) #not to overload your computer
# registerDoParallel(cl)

##################################################################################

file_exptInfo <- here::here("data", "referenceData/sample_info.txt")
geneCdsFile <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_version_s10-m04-r03_CDS_Unique.bed"
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
orgDb <- org.Anidulans.eg.db


TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")
other_dataPath <- here::here("data", "other_data")

matrixType <- "deeptools"
matrixDim = c(200, 200, 100, 10)

geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::select(-score) %>% 
  dplyr::mutate(length = end - start)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID"))

print(TF_dataPath)
print(file_exptInfo)
##################################################################################


polIIsampleFile <- paste(polII_dataPath, "/", "sample_polII.list", sep = "")

polIISampleList <- fread(file = polIIsampleFile, sep = "\t", header = F,
                         stringsAsFactors = F, col.names = c("id"), data.table = F)


polII_info <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = polIISampleList$id,
  dataPath = polII_dataPath,
  matrixSource = matrixType)

i <- 1

## process all polII expression matrix
foreach(i = 1:nrow(polII_info),
        .packages = c("chipmine"),
        .verbose = TRUE) %dopar% {
          
          polIIDf <- chipmine::preProcess_polII_expression(
            expMat = polII_info$polIIExpMat[i],
            title = polII_info$sampleId[i],
            expFraction = 10,
            polIIExpFile = polII_info$polIIExpFile[i])
          
          # mat5Kb <- gsub(pattern = ".tab.gz", replacement = "_5kb.tab.gz", x = polII_info$matFile[i])
          # ## 5kb - 2kb - 1kb matrix
          # bwMat <- chipmine::bigwig_profile_matrix(
          #   bwFile = polII_info$bwFile[i],
          #   bedFile = file_genes,
          #   signalName = polII_info$sampleId[i],
          #   genes = geneSet$gene,
          #   readLocal = FALSE,
          #   storeLocal = TRUE,
          #   localPath = mat5Kb,
          #   extend = c(5000, 1000),
          #   target_ratio = 0.25)
          
          polII_info$bwFile[i]
        }


##################################################################################
## process all TF macs2 results to merge everything
tfSampleFile <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")
# tfSampleFile <- paste(TF_dataPath, "/", "sample_tfs.list", sep = "")

tfSampleList <- readr::read_tsv(file = tfSampleFile, col_names = c("id"),  comment = "#") %>%
  as.data.frame()


tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSampleList$id,
  dataPath = TF_dataPath,
  matrixSource = matrixType)



i <- 1
# foreach(i = 1:nrow(tfInfo),
#         .packages = c("chipmine"),
#         .verbose = TRUE) %do% {

for(i in 1:nrow(tfInfo)){

  peakAn <- narrowPeak_annotate(
    peakFile = tfInfo$narrowpeakFile[i],
    txdb = txDb,
    includeFractionCut = 0.7,
    bindingInGene = FALSE,
    promoterLength = 3000,
    insideSkewToEndCut = 0.7,
    reportPseudo = FALSE,
    output = tfInfo$narrowpeakAnno[i])
  
  
  tfDf <- gene_level_peak_annotation(
    sampleId = tfInfo$sampleId[i],
    peakAnnotation = tfInfo$narrowpeakAnno[i],
    genesDf = geneSet,
    peakFile = tfInfo$narrowpeakFile[i],
    bwFile = tfInfo$bwFile[i],
    outFile = tfInfo$peakTargetFile[i],
    bindingInGene = FALSE)
  
  # ## preprocess data
  # tfDf <- chipmine::preProcess_macs2_results(
  #   title = tfInfo$sampleId[i],
  #   peakAnnoFile = tfInfo$narrowpeakAnno[i],
  #   cdsFile = geneCdsFile,
  #   peakFile = tfInfo$narrowpeakFile[i],
  #   bwFile = tfInfo$bwFile[i],
  #   outFile = tfInfo$tfPeakFile[i],
  #   bindingInGene = FALSE)
  # 
  # ## plot for all genes
  # outPath <- paste(TF_dataPath, "/", tfInfo$sampleId[i], sep = "")
  # outPrefix_all <- paste0(outPath, "/", tfInfo$sampleId[i], "_allGenes", collapse = "")
  # mat <- chipmine::import_profile_from_file(
  #   file = tfInfo$matFile[i],
  #   source = matrixType,
  #   signalName = tfInfo$sampleId[i],
  #   selectGenes = geneSet$gene
  # )
  # 
  # profCol <- colorRamp2(quantile(mat, c(0.50, 0.995), na.rm = T), c("white", "red"))
  # ht1 <- chipmine::profile_heatmap(
  #   profileMat = mat,
  #   signalName = tfInfo$profileName[i],
  #   profileColor = profCol,
  #   columnTitle = tfInfo$sampleId[i],
  #   geneGroups = tfInfo$clusterFile[i]
  # )
  # 
  # anGl <- gene_length_heatmap_annotation(bedFile = file_genes, genes = geneSet$gene)
  # 
  # htlist_allGenes <- anGl$an + ht1$rowAnno + ht1$heatmap
  # 
  # all_title <- paste0(tfInfo$sampleId[i], ": all genes", collapse = "")
  # 
  # pdf(file = paste0(outPrefix_all, ".pdf", collapse = ""), width = 8, height = 12)
  # draw(htlist_allGenes,
  #      main_heatmap = tfInfo$profileName[i],
  #      column_title = all_title,
  #      column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  #      row_sub_title_side = "left",
  #      gap = unit(5, "mm"),
  #      padding = unit(rep(0.5, times = 4), "cm")
  # )
  # 
  # 
  # row_annotation_axis(an = "gene_length",
  #                     at = c(0, 2000, 4000),
  #                     labels = c("0kb", "2kb", ">4kb"),
  #                     slice = length(unique(ht1$cluster$cluster)))
  # 
  # 
  # dev.off()
  # 
  
  
  # mat5Kb <- gsub(pattern = ".tab.gz", replacement = "_5kb.tab.gz", x = tfInfo$matFile[i])
  # 
  # ## 5kb - 2kb - 1kb matrix
  # bwMat <- chipmine::bigwig_profile_matrix(
  #   bwFile = tfInfo$bwFile[i],
  #   bedFile = file_genes,
  #   signalName = tfInfo$sampleId[i],
  #   genes = geneSet$gene,
  #   readLocal = FALSE,
  #   storeLocal = TRUE,
  #   localPath = mat5Kb,
  #   extend = c(5000, 1000),
  #   target_ratio = 0.25)
  
  # tfInfo$bwFile[i]
  
}



##################################################################################

# stopCluster(cl)


