# library(ballgown)
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(FactoMineR))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))

## This script uses tximport to extract the FPKM matrix from the stringTie output
## It also perform PCA analysis using FPKM values
rm(list = ls())

analysisName <- "count_data"

file_sampleInfo <- here::here("data", "reference_data", "polII_sample_info.txt")
readLength <- as.numeric(readr::read_file(file = here::here("data", "read_length.config")))

outDir <- here("analysis", "02_QC_polII")
outPrefix <- paste(outDir, "/", analysisName, sep = "")
orgDb <- org.Anidulans.FGSCA4.eg.db
col_geneId <- "GID"

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}


###########################################################################
exptInfo <- suppressMessages(readr::read_tsv(file = file_sampleInfo)) %>% 
  as.data.frame()

rownames(exptInfo) <- exptInfo$sampleId

## set the factor levels.
## control levels should be first
# exptInfo$gt <- forcats::as_factor(exptInfo$gt)
# exptInfo$treatment <- forcats::as_factor(exptInfo$treatment)
exptInfo$condition <- forcats::as_factor(exptInfo$condition)

## Add gene symbol for each Ensembl ID
geneInfo <- AnnotationDbi::select(x = orgDb,
                                  keys = keys(x = orgDb, keytype = col_geneId),
                                  columns = c("GENE_NAME", "DESCRIPTION"),
                                  keytype = col_geneId) %>% 
  dplyr::rename(geneId = !!col_geneId)

design <- ~ condition

###########################################################################
## import counts data: either by tximport or as raw count matrix

# ## import the counts data using tximport and run DESeq2
# path_stringtie <- here::here("data", "stringTie")
# filesStringtie <- paste(path_stringtie, "/stringTie_", exptInfo$sampleId, "/t_data.ctab", sep = "")
# names(filesStringtie) <- exptInfo$sampleId
# 
# tmp <- data.table::fread(file = filesStringtie[1], sep = "\t", header = T, stringsAsFactors = F)
# tx2gene <- tmp[, c("t_name", "gene_id")]
# 
# cat("Importing stringTie data with tximport. Read length = ", readLength)
# 
# txi <- tximport(files = filesStringtie, type = "stringtie",
#                 tx2gene = tx2gene, readLength = readLength)
# 
# ddsTxi <- DESeqDataSetFromTximport(txi = txi, colData = exptInfo, design = design)
# # assay(ddsTxi)
# # colData(ddsTxi)
# # rowData(ddsTxi)
# 
# ## Run DESeq2
# dds <- DESeq(ddsTxi)

## import raw counts data and run DESeq2
file_rawCounts <- here::here("data", "polII_data", "polII_raw_counts.tab")

countsDf <- suppressMessages(readr::read_tsv(file = file_rawCounts, col_names = T)) %>%
  as.data.frame()
rownames(countsDf) <- countsDf$geneId
countsDf$geneId <- NULL


if(all(rownames(exptInfo) %in% colnames(countsDf))){
  countsDf <- countsDf[, rownames(exptInfo)]
} else{
  stop("Column names in count matrix does not match with row names in experiment data")
}

## select only those sample rows which are part of current comparison
# exptInfo <- droplevels(subset(exptInfo, condition %in% compare))
# countsDf <- countsDf[, rownames(exptInfo)]

## run DESeq2 and extract the processed data
ddsCount <- DESeqDataSetFromMatrix(countData = countsDf, colData = exptInfo, design = design)

## Run DESeq2
dds <- DESeq(ddsCount)

###########################################################################
## raw counts
# rawCounts <- tibble::rownames_to_column(as.data.frame(counts(dds, normalized = FALSE)), var = "geneId")
# 
# readr::write_tsv(x = rawCounts, path = paste0(c(outPrefix,".rawCounts.tab"), collapse = ""))

# ## FPKM
# fpkmCounts <- tibble::rownames_to_column(as.data.frame(fpkm(dds)), var = "geneId");
# 
# if(all(fpkmCounts$geneId %in% geneInfo$geneId)){
#   fpkmCounts <- dplyr::left_join(x = fpkmCounts, y = geneInfo, by = c("geneId" = "geneId")) %>% 
#     dplyr::select(geneId, exptInfo$sampleId, everything()) %>% 
#     dplyr::filter(!is.na(geneId))
#   
#   readr::write_tsv(x = fpkmCounts, path = paste0(c(outPrefix,".FPKM.tab"), collapse = ""))
#   
# } else{
#   warning(
#     "Missing geneIds in fpkmCounts from org.db:\n",
#     paste(head(fpkmCounts$geneId[which(!fpkmCounts$geneId %in% geneInfo$geneId)], 10), collapse = ", "),
#     ", ..."
#   )
# }



## normalized counts matrix
normCounts <- tibble::rownames_to_column(as.data.frame(counts(dds, normalized = TRUE)), var = "geneId")
readr::write_tsv(x = normCounts, path = paste(outPrefix,".normCounts.tab", sep = ""))

## r-log normalized counts
# rld <- rlog(dds, blind = FALSE)
rld <- vst(object = dds)
rldCount <- rownames_to_column(as.data.frame(assay(rld)), var = "geneId")

readr::write_tsv(x = rldCount, path = paste(outPrefix,".rlogCounts.tab", sep = ""))

plotPCA(rld, intgroup = "condition", ntop = 4000)

pcaData <- plotPCA(rld, intgroup= "condition", returnData = TRUE, ntop = 4000)
percentVar <- sprintf("%.2f", 100 * attr(pcaData, "percentVar"))

# pcaData$treatment <- forcats::as_factor(pcaData$treatment)
# pcaData$genotype <- forcats::as_factor(pcaData$genotype)
pcaData$condition <- forcats::as_factor(pcaData$condition)

pltTitle <- "Principal Component Analysis"
# pointCol <- base::structure(RColorBrewer::brewer.pal(n = length(unique(pcaData$condition)), name = "Set1"),
#                             names = levels(pcaData$condition))


pcaPlot <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(mapping = aes(color = condition), shape = 16, size=4) +
  geom_text_repel(mapping = aes(label = name), size = 3, point.padding = 0.5) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  # scale_color_manual(values = pointCol) +
  ggtitle(pltTitle) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position = "none",
        legend.text = element_text(size = 13),
        legend.title = element_text(face = "bold", size = 15)
  )


# png(filename = paste(outPrefix, ".PCA.png", sep = ""), width = 10000, height = 10000, res = 700)
pdf(file = paste(outPrefix, ".PCA.pdf", sep = ""), width = 15, height = 15)
pcaPlot
dev.off()


###########################################################################
## sample distance matrix

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)


plot_dist <- ComplexHeatmap::Heatmap(
  matrix = sampleDistMatrix,
  col = colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255),
  column_title = "Distance matrix of normalized read counts",
  column_title_gp = gpar(fontface = "bold", fontsize = 14),
  row_names_max_width = unit(10, "cm"),
  column_names_max_height = unit(10, "cm"),
  row_dend_width = unit(3, "cm"),
  column_dend_height = unit(3, "cm"),
  heatmap_legend_param = list(title = "Distance", title_gp = gpar(fontsize = 12),
                              title_position = "topcenter")
)


png(filename = paste(outPrefix, ".distance_heatmap.png", sep = ""),
    width = 10000, height = 10000, res = 350)
draw(plot_dist,
     padding = unit(rep(0.5, 4), "cm")
)
dev.off()


###########################################################################
## read the FPKM data from stringTie output using ballgown package
## no clearity on how ballgown calculate gene level FPKM scores using transcript score.
## so use stringTie -> tximport way
## tximport calculate gene level FPKM/TPM scores by using weighted average over all transcripts
## for a gene. See tximport code for the details
# filesStringtie <- paste(path_stringtie, "/stringTie_", exptInfo$sampleId, sep = "")
# 
# bg <- ballgown(samples = filesStringtie, meas='all')
# 
# geneExpression <- data.frame(gexpr(bg)) %>%
#   rownames_to_column(var = "geneId")
# 
# 
# ## rename the sample names
# renameCols <- base::structure(gsub(pattern = "FPKM.stringTie_", replacement = "", x = names(geneExpression)),
#                               names = names(geneExpression))
# 
# 
# renameCols[colnames(geneExpression)]
# 
# colnames(geneExpression) <- renameCols[colnames(geneExpression)]
# 
# 
# fwrite(x = geneExpression, file = paste(outPrefix, "FPKM_ballgown.tab", sep = "_"),
#        sep = "\t", col.names = T, row.names = F, quote = F)
# 
# 



###########################################################################
## PCA based on rld counts

normCountMat <- as.matrix(rldCount[, c(exptInfo$sampleId), drop = FALSE])
rownames(normCountMat) <- rldCount$geneId

## remove low count rows
keep <- rowSums(normCountMat > 1) >= 2
normCountMatFiltered <- normCountMat[keep, ]

## transform the data such that the genes are columns and each sample is a row
## also append the additional information for each sample using left_join()
df2 <- as.data.frame(t(normCountMatFiltered)) %>% 
  tibble::rownames_to_column(var = "sampleId") %>% 
  dplyr::left_join(y = exptInfo, by = c("sampleId" = "sampleId")) %>% 
  dplyr::select(!!!colnames(exptInfo), dplyr::everything())

row.names(df2) <- df2$sampleId

res.pca <- PCA(df2, graph = FALSE, scale.unit = TRUE,
               quali.sup = 1:ncol(exptInfo), ncp = 10)

eig.val <- get_eigenvalue(res.pca)

## scree plot: variance by PC
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

## Graph of individuals
ind <- get_pca_ind(res.pca)

fviz_pca_ind(res.pca,
             # col.ind = exprData$treatment,
             fill.ind = df2$condition,
             pointshape = 21,
             repel = TRUE,
             mean.point = FALSE,
             legend.title = "Study",
             pointsize = 3
)


## prepare the plot dataframe for ggplot
plotData <- as.data.frame(ind$coord) %>%
  tibble::rownames_to_column(var = "sampleId") %>%
  dplyr::left_join(y = exptInfo, by = c("sampleId" = "sampleId"))

pairs(x = plotData[, 2:6],
      pch = 19,  cex = 1,
      col = plotData$condition,
      lower.panel=NULL)


## set the factor levels
plotData$condition <- forcats::as_factor(plotData$condition)
plotData$genotype <- forcats::as_factor(plotData$genotype)
plotData$treatment <- forcats::as_factor(plotData$treatment)

pointCol <- base::structure(
  RColorBrewer::brewer.pal(n = length(unique(plotData$condition)), name = "Set1"),
  names = levels(plotData$condition))

pltTitle <- "Principal Component Analysis"


## decide which PCs to use for plotting
pcToPlot <- c(1, 2)
pcCols <- grep(pattern = "Dim.", x = colnames(plotData), value = T)[pcToPlot]
fillColumn <- "condition"
# shapeColumn <- "treatment"

pcaPlot <- ggplot(data = plotData,
                  mapping = aes(x = !!sym(pcCols[1]), y = !!sym(pcCols[2]), label = sampleId)) +
  geom_point(mapping = aes(color = !!sym(fillColumn)),
             size = 4, stroke = 2) +
  # scale_shape_manual(values = c(1, 15, 17)) +
  # guides(fill=guide_legend(override.aes=list(shape=21))) +
  scale_color_manual(values = pointCol) +
  geom_text_repel(size = 3, point.padding = 0.5) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) + 
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
  xlab( paste("PC",pcToPlot[1]," (", sprintf("%.2f", eig.val[pcToPlot[1], "variance.percent"]), "%)", sep = "") ) +
  ylab( paste("PC",pcToPlot[2]," (", sprintf("%.2f", eig.val[pcToPlot[2], "variance.percent"]), "%)", sep = "") ) +
  ggtitle(pltTitle) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold"))



# pdf(file = paste(outPrefix, ".rld_PCA.pdf", sep = ""), width = 10, height = 10)
png(filename = paste(outPrefix, ".rld_PCA.png", sep = ""), width = 4000, height = 3000, res = 350)
print(pcaPlot)
dev.off()

#############################################################################
## correlation scatter plot
pt <- ggpairs(
  data = dplyr::select(rldCount, -geneId),
              upper = list(continuous = wrap("points", size = 0.1)),
              lower = list(continuous = wrap("cor", size = 10)),
              diag = list(continuous = "densityDiag"),
  title = "scatter plot of rlog transformed normalized read counts") +
  theme_bw() +
  theme(
    strip.text.y = element_text(size = 18, angle = 0, hjust = 0, face = "bold"),
    strip.text.x = element_text(size = 18, angle = 90, hjust = 0, face = "bold")
  )

png(filename = paste(outPrefix, ".scatter_matrix.png", sep = ""),
    width = 10000, height = 10000, res = 300)

pt
dev.off()

# ###########################################################################
# ## calculate mean FPKM
# 
# exptInfo <- read.table(file = sampleInfoFile, header = T, sep = "\t", row.names = "sampleId")
# 
# ## IF NEEDED: modify the row names and column names
# # colnames(rawData) <- sub("_WT", "_", colnames(rawData))
# # rownames(exptInfo) <- sub("_WT", "_", rownames(exptInfo))
# 
# ## select only those sample rows which are part of current comparison
# designInfo <- exptInfo
# # designInfo <- droplevels(subset(exptInfo, condition %in% compare))
# 
# if(! all( rownames(designInfo) %in% colnames(geneExpression) ) ){
#   stop("Column names in FPKM matrix does not match with row names in experiment data")
# }
# 
# 
# ## get the sample names for each condition under comparison
# grp1 <- sapply(rownames(designInfo[designInfo$condition %in% compare[1], ]), FUN = as.name, USE.NAMES = F, simplify = T)
# name1 <- paste(compare[1], "_meanFPKM", sep = "")
# grp1Len <- length(grp1)
# 
# grp2 <- sapply(rownames(designInfo[designInfo$condition %in% compare[2], ]), FUN = as.name, USE.NAMES = F, simplify = T)
# name2 <- paste(compare[2], "_meanFPKM", sep = "")
# grp2Len <- length(grp2)
# 
# ## can add multiple groups
# 
# ## calculate the average FPKM
# fpkmData <- geneExpression %>% 
#   dplyr::select(geneId, !!!c(grp1, grp2)) %>%
#   rowwise() %>%
#   mutate(!!name1 := sum(!!!grp1) / !!grp1Len,
#          !!name2 := sum(!!!grp2) / !!grp2Len)
# 
# 
# 
# ## write the data to file
# fwrite(x = fpkmData, file = "FPKM_matrix.tab", sep = "\t", quote = F, col.names = T)
# 
# 


