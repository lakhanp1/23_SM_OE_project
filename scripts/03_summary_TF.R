library(chipmine)
library(org.Anidulans.eg.db)
library(TxDb.Anidulans.AspGD.GFF)
library(here)
library(ggridges)
library(ggbeeswarm)
library(ggpubr)
library(ggrepel)

## 1) peak enrichment distribution
## 2) peak p-value distribution
## 3) peak annottion pie chart
## 4) combined matrix of peak enrichment distribution plots
## 5) combined matrix of peak p-value distribution plots

rm(list = ls())

##################################################################################
analysisName <- "TF_ChIP_summary"
outDir <- here::here("analysis", "03_QC_TF")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_exptInfo <- here::here("data", "referenceData/sample_info.txt")

file_genes <- here::here("data", "referenceData/AN_genes_for_polII.bed")
orgDb <- org.Anidulans.eg.db
txDb <- TxDb.Anidulans.AspGD.GFF

TF_dataPath <- here::here("data", "TF_data")

file_tf_macs2 <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")
file_tf <- paste(TF_dataPath, "/", "sample_tf.list", sep = "")

matrixType <- "2kb_summit"
up <- 2000
down <- 2000
body <- 0
bin <- 10
matrixDim = c(c(up, body, down)/bin, bin)

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

theme_scatter <- theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    panel.grid = element_blank()
  )

allPlotData <- NULL

pdf(file = paste(outPrefix, ".macs2.pdf", sep = ""), width = 10, height = 12, onefile = TRUE)

i <- 98

for (i in 1:nrow(tfInfo)) {
  
  print(tfInfo$sampleId[i])
  
  ## annotate peaks and prepare gene level annotation file
  peakType <- dplyr::case_when(
    tfInfo$peakType[i] == "narrow" ~ "narrowPeak",
    tfInfo$peakType[i] == "broad" ~ "broadPeak"
  )
  
  peaksGr <- rtracklayer::import(con = tfInfo$peakFile[i], format = peakType)
  
  if (file.exists(tfInfo$peakAnno[i])) {
    
    backboneGene <- tfInfo$SM_TF[i]
    
    ## few TFs are mapped to multiple SM clusters. So preparing the list of SM tf data
    smTfInfo <- suppressMessages(
      AnnotationDbi::select(
        x = orgDb,
        keys = na.omit(backboneGene),
        columns = c("SM_GENE", "SM_CLUSTER", "SM_ID"),
        keytype = "GID")) %>% 
      dplyr::filter(!is.na(SM_ID))
    
    clusterGenes <- suppressMessages(
      AnnotationDbi::select(x = orgDb,
                            keys = smTfInfo$SM_ID,
                            columns = c("GID", "SM_GENE"),
                            keytype = "SM_ID")
    )
    
    genesToMark <- list(
      SM_cluster = clusterGenes$SM_GENE,
      other_SM_genes = setdiff(smGenes$SM_GENE, clusterGenes$SM_GENE)
    )
    
    markingPreference <- tibble(geneType = c(names(genesToMark), "peaks"),
                                drawOrder = 1:(length(genesToMark)+1))
    
    peakMarkDf <- purrr::map_dfr(
      .x = genesToMark,
      .f = function(x){data.frame(geneId = x, stringsAsFactors = F)},
      .id = "geneType")
    
    ## import peak data with annotation
    peakAnno <- import_peak_annotation(sampleId = tfInfo$sampleId[i],
                                       peakAnnoFile = tfInfo$peakAnno[i],
                                       renameColumn = FALSE)
    
    
    ## add target gene type information if points need to be colored
    peakAnno <- dplyr::left_join(x = peakAnno, y = peakMarkDf, by = c("geneId" = "geneId")) %>% 
      tidyr::replace_na(replace = list(geneType = "peaks")) %>% 
      dplyr::left_join(y = markingPreference, by = c("geneType" = "geneType"))
    
    ## for each peak, select one gene. preference is decided by order of names in genesToMark list
    plotData <- dplyr::group_by(peakAnno, peakId) %>% 
      dplyr::arrange(drawOrder) %>% 
      dplyr::slice(1L) %>%
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(sampleId = tfInfo$sampleId[i])
    
    plotData$geneType <- factor(x = plotData$geneType,
                                levels = markingPreference$geneType,
                                ordered = TRUE)
    
    
    ## summary table
    summaryTable <- purrr::map_dfr(
      .x = structure(c("peakEnrichment", "peakPval", "peakQval"),
                     names = c("peakEnrichment", "peakPval", "peakQval")),
      .f = function(x){
        as.list(summary(plotData[[x]]))
      },
      .id = "metric")
    
    gg_stable <- ggtexttable(summaryTable, rows = NULL, 
                             theme = ttheme("mOrange"))
    
    
    ## set outliers to 99.5 quantile
    if(nrow(plotData) > 20){
      plotData$peakEnrichment <- pmin(plotData$peakEnrichment, quantile(plotData$peakEnrichment, 0.995))
      plotData$peakPval <- pmin(plotData$peakPval, quantile(plotData$peakPval, 0.99))
    }
    
    allPlotData <- dplyr::bind_rows(allPlotData, plotData)
    
    quantile(plotData$peakEnrichment,
             c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
    
    quantile(plotData$peakPval,
             c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
    
    pointColor <- structure(c("red", "blue", "grey"), names = markingPreference$geneType)
    pointAlpha <- structure(c(1, 0.5, 0.5), names = markingPreference$geneType)
    
    
    ## macs2 fold-enrichment density scatter plot
    gg_dot_enrichment <- ggplot(
      data = plotData,
      mapping = aes(x = tfInfo$sampleId[i], y = peakEnrichment)) +
      geom_quasirandom(mapping = aes(color = geneType, alpha = geneType)) +
      geom_boxplot(width=0.1, fill = NA, outlier.colour = NA, color = alpha("black", 0.7)) +
      scale_color_manual(name = "Peak annotations",
                         values = pointColor) +
      scale_alpha_manual(values = pointAlpha) +
      labs(title = "macs2 fold enrichment distribution",
           y = "peak enrichment") +
      guides(alpha = FALSE,
             color = guide_legend(override.aes = list(size = 4))) +
      theme_scatter
    
    
    
    ## macs2 p-value density scatter plot
    gg_dot_pval <- ggplot(
      data = plotData,
      mapping = aes(x = tfInfo$sampleId[i], y = peakPval)) +
      geom_quasirandom(mapping = aes(color = geneType, alpha = geneType)) +
      geom_boxplot(width=0.1, fill = NA, outlier.colour = NA, color = alpha("black", 0.7)) +
      scale_color_manual(name = "Peak annotations",
                         values = pointColor) +
      scale_alpha_manual(values = pointAlpha) +
      labs(title = "macs2 p-value distribution",
           y = "-log10(p-value)") +
      guides(alpha = FALSE,
             color = guide_legend(override.aes = list(size = 4))) +
      theme_scatter
    

    ## peak annotation pie chart
    peakAnSummary <- dplyr::group_by(peakAnno, peakType) %>% 
      dplyr::summarise(count = n()) %>% 
      dplyr::mutate(label = round(count / sum(count), digits = 4)) %>% 
      dplyr::mutate(
        peakType = factor(peakType,
                          levels = c("upstream", "promoter", "include_tx", "5UTR", "tx_start",
                                     "EXON", "INTRON", "tx_end", "3UTR", "intergenic"))
      )
    
    peakTypeCol <- c("upstream" = "#a6cee3", "promoter" = "#1f78b4", "include_tx" = "#b2df8a",
                     "5UTR" = "#fb9a99", "tx_start" = "#e31a1c", "EXON" = "#ff7f00",
                     "INTRON" = "#fdbf6f", "tx_end" = "#cab2d6", "3UTR" = "#6a3d9a",
                     "intergenic" = "#b15928")
    
    
    gg_pie_peakAn <- ggplot(data = peakAnSummary,
                            mapping = aes(x = 1, y = count, fill = peakType)) +
      geom_bar(color = "black", stat = "identity") +
      coord_polar(theta = "y", start = 0, direction = -1) +
      geom_label_repel(mapping = aes(x = 1.4, label = scales::percent(label)),
                       position = position_stack(vjust = 0.5),
                       size = 5, show.legend = FALSE) +
      scale_fill_manual(
        name = "Peak annotations",
        values = peakTypeCol
      ) +
      labs(title = "Peak annotation distribution") +
      theme_bw() +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title=element_text(size=14, face="bold")
      )
    
    
    ## combine plots + table to make a summary figure
    summaryFig <- ggarrange(
      gg_stable,
      ggarrange(gg_dot_enrichment, gg_dot_pval,
                nrow = 1, ncol = 2, legend = "bottom", common.legend = TRUE),
      gg_pie_peakAn,
      nrow = 3,
      legend = "right",
      heights = c(1, 4, 4)
    )
    
    
    summaryFig <- annotate_figure(
      p = summaryFig,
      top = text_grob(label = paste(tfInfo$sampleId[i], "ChIPseq summary\n#peaks=",length(peaksGr)),
                      size = 16, face = "bold")
    )
    
  } else{
    
    summaryFig <- annotate_figure(
      p = ggplot() + geom_text(mapping = aes(x = 0.5, y = 0.5), label = "No data", size = 30) +
        theme_void(),
      top = text_grob(label = paste(tfInfo$sampleId[i], "ChIPseq summary\n#peaks =",length(peaksGr)),
                      size = 16, face = "bold")
    )
    
  }
  
  plot(summaryFig)
}

dev.off()

##################################################################################
## combined summary plot matrix

theme_plot_matrix <- theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    title = element_text(hjust = 0.5, size = 20),
    strip.text = element_text(size = 8, hjust = 0),
    strip.background = element_rect(fill = "white")
  )


## combined summary plot matrix: peak enrichment
gg_all_enrichment <- ggplot(
  data = allPlotData,
  mapping = aes(x = sampleId[i], y = peakEnrichment)) +
  geom_quasirandom(color = "#bf9a2d") +
  geom_boxplot(width=0.1, fill = NA, outlier.colour = NA, color = alpha("black", 1)) +
  geom_hline(yintercept = 3, color = "blue") +
  labs(title = "masc2 fold enrichment distribution") +
  facet_wrap(facets = ~ sampleId, ncol = 14, scales = "free") +
  theme_plot_matrix

# pdf(file = paste(outPrefix, ".macs2_enrichment.pdf", sep = ""), width = 16, height = 16, onefile = TRUE)
png(filename = paste(outPrefix, ".macs2_enrichment.png", sep = ""), width = 15000, height = 10000, res = 350)
gg_all_enrichment
dev.off()

## combined summary plot matrix: peak p-value
gg_all_pval <- ggplot(
  data = allPlotData,
  mapping = aes(x = sampleId[i], y = peakPval)) +
  geom_quasirandom(color = "#567a0f") +
  geom_boxplot(width=0.1, fill = NA, outlier.colour = NA, color = alpha("black", 1)) +
  geom_hline(yintercept = 20, color = "red") +
  labs(title = "masc2 p-value distribution") +
  facet_wrap(facets = ~ sampleId, ncol = 14, scales = "free") +
  theme_plot_matrix

# pdf(file = paste(outPrefix, ".macs2_enrichment.pdf", sep = ""), width = 16, height = 16, onefile = TRUE)
png(filename = paste(outPrefix, ".macs2_pval.png", sep = ""), width = 15000, height = 10000, res = 350)
gg_all_pval
dev.off()


##################################################################################




##################################################################################
# ## profile plot
# 
# outDir <- dirname(tfInfo$matFile[i])
# outPrefix <- paste(outDir, "/", tfInfo$sampleId[i], ".raw_peaks_summary", sep = "")
# 
# ## control samples to plot alongside TF ChIP sample
# tfControls <- c("AN10300_sCopy_OE_16h_input_FLAG_ChIPMix55_1",
#                 "AN10300_sCopy_OE_16h_input_FLAG_ChIPMix55_2",
#                 "AN0148_sCopy_OE_16h_HA_ChIPMix46_1",
#                 "AN2025_sCopy_OE_16h_HA_ChIPMix36_1")
# 
# 
# ## create profile matrix of 2kb region around peak summit for control samples
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
#   ctrlSampleInfo <- get_sample_information(
#     exptInfoFile = file_exptInfo,
#     samples = tfControls,
#     dataPath = TF_dataPath,
#     profileMatrixSuffix = matrixType)
#   
#   for (ctrlIdx in 1:nrow(ctrlSampleInfo)) {
#     profileMat <- bigwig_profile_matrix(bwFile = ctrlSampleInfo$bwFile[ctrlIdx],
#                                         regions = peakSummitGr,
#                                         signalName = ctrlSampleInfo$profileName[ctrlIdx],
#                                         extend = c(up, down),
#                                         targetName = "summit",
#                                         storeLocal = T,
#                                         localPath = ctrlSampleInfo$matFile[ctrlIdx])
#   }
#   
#   
#   
#   tempSampleInfo <- dplyr::bind_rows(tfInfo[i, ], ctrlSampleInfo) %>% 
#     dplyr::distinct()
#   
#   tfProfileMat <- import_profile_from_file(
#     file = tfInfo$matFile[i],
#     signalName = tfInfo$profileName[i],
#     selectGenes = peakSummitGr$name,
#     up = matrixDim[1], target = matrixDim[2], down = matrixDim[3],
#     targetType = "point", targetName = "summit" 
#   )
#   
#   quantile(tfProfileMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
#   # tfMeanColor <- colorRamp2(quantile(tfProfileMat, c(0.50, 0.995), na.rm = T), c("white", "red"))
#   tfColorList <- sapply(
#     X = tempSampleInfo$sampleId,
#     FUN = function(x){
#       return(colorRamp2(breaks = quantile(tfProfileMat, c(0.50, 0.995), na.rm = T),
#                         colors = c("white", "red")))
#     }
#   )
#   
#   
#   
#   profilePlots <- multi_profile_plots(exptInfo = tempSampleInfo,
#                                       genesToPlot = peakSummitGr$name,
#                                       profileColors = tfColorList,
#                                       targetType = "point", targetName = "summit",
#                                       clusters = NULL,
#                                       showAnnotation = FALSE,
#                                       matBins = matrixDim,
#                                       column_title_gp = gpar(fontsize = 12))
#   
#   
#   rowOrd_peaks <- order(peakSummitGr$signalValue, decreasing = TRUE)
#   
#   # pdf(file = paste(outPrefix, "_profiles.pdf", sep = ""), width = 18, height = 13)
#   png(file = paste(outPrefix, ".profiles.png", sep = ""), width = 3500, height = 2500, res = 250)
#   
#   ht <- draw(
#     profilePlots$heatmapList,
#     main_heatmap = tempSampleInfo$profileName[1],
#     row_order = rowOrd_peaks,
#     column_title = paste(tfInfo$sampleId[i], "binding signal around 2kb region of macs2 peak summit"),
#     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
#     row_sub_title_side = "left",
#     heatmap_legend_side = "bottom",
#     gap = unit(7, "mm"),
#     padding = unit(rep(0.5, times = 4), "cm")
#   )
#   
#   dev.off()
#   
# }
# 
# 










