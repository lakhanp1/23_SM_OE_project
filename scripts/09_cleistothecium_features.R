library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(TxDb.Anidulans.FGSCA4.AspGD.GFF)
library(FactoMineR)
library(factoextra)

## data mining to extract imp genes with role in cleistothecium formation

rm(list = ls())

##################################################################################

analysisName <- "cleistothecium"
outDir <- here::here("analysis", "05_cleistothecium")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_classes <- paste(outDir, "/ascocarp_classes.txt", sep = "")
file_exptInfo <- here::here("data", "referenceData/sample_info.txt")
file_sexAsex <- here::here("data", "referenceData/sex_asex_genes.txt")

file_genes <- here::here("data", "referenceData/AN_genes_for_polII.bed")
orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")

file_polIISamples <- paste(polII_dataPath, "/", "sample_polII.list", sep = "")

geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "geneId", "score", "strand")) %>%
  dplyr::select(geneId)


# geneDesc <- select(x = orgDb, keys = geneSet$geneId,
#                    columns = c("GENE_NAME", "DESCRIPTION"), keytype = "GID")
# 
# geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID"))

##################################################################################

sexasexGenes <- suppressMessages(readr::read_tsv(file = file_sexAsex))

polIISampleList <- readr::read_tsv(file = file_polIISamples, col_names = c("sampleId"),  comment = "#")

polII_info <- get_sample_information(exptInfoFile = file_exptInfo,
                                     samples = polIISampleList$sampleId,
                                     dataPath = polII_dataPath,
                                     profileMatrixSuffix = "normalizedmatrix")

polII_ids <- polII_info$sampleId[which(polII_info$IP_tag == "polII")]

polIICols <- list(exp = structure(polII_ids, names = polII_ids)) %>% 
  append(
    sapply(X = c("is_expressed", "transformed"),
           FUN = function(x){ structure(paste(x, ".", polII_ids, sep = ""), names = polII_ids) },
           simplify = F, USE.NAMES = T)
  )


sampleClasses <- suppressMessages(readr::read_tsv(file = file_classes))
sampleClasses$ascocarp <- forcats::as_factor(sampleClasses$ascocarp)

polIIMat <- get_polII_expressions(genesDf = geneSet, exptInfo = polII_info)


quantRankFormula <- sapply(
  X = polII_ids,
  FUN = function(x){
    rhs <- paste("transformed", ".", x, sep = "")
    form <- quos(!!rhs := dplyr::percent_rank(!!sym(x)))
    return(form)
  },
  simplify = F, USE.NAMES = F
) %>% 
  unlist()

polIIMat <- dplyr::mutate(.data = polIIMat, !!!quantRankFormula)


##################################################################################

fpkmMat <- as.matrix(polIIMat[, polIICols$transformed, drop = FALSE])
rownames(fpkmMat) <- polIIMat$geneId
colnames(fpkmMat) <- stringr::str_replace(string = colnames(fpkmMat),
                                          pattern = "transformed.", replacement = "")

fpkmMat <- fpkmMat[sexasexGenes$geneId, ]

df2 <- as.data.frame(t(fpkmMat)) %>% 
  tibble::rownames_to_column(var = "sampleId") %>% 
  dplyr::left_join(y = sampleClasses, by = "sampleId") %>% 
  dplyr::select(!!!colnames(sampleClasses), dplyr::everything())

row.names(df2) <- df2$sampleId

res.pca <- PCA(df2, graph = FALSE, scale.unit = TRUE,
               quali.sup = 1:ncol(sampleClasses), ncp = 10)

eig.val <- get_eigenvalue(res.pca)

## scree plot: variance by PC
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

## Graph of individuals
ind <- get_pca_ind(res.pca)

fviz_pca_ind(res.pca,
             # col.ind = exprData$treatment,
             fill.ind = df2$ascocarp,
             geom = "point",
             pointshape = 21,
             repel = TRUE,
             mean.point = FALSE,
             legend.title = "cleistothecium",
             pointsize = 3
)


## prepare the plot dataframe for ggplot
plotData <- as.data.frame(ind$coord) %>%
  tibble::rownames_to_column(var = "sampleId") %>%
  dplyr::left_join(y = sampleClasses, by = c("sampleId" = "sampleId"))


pairs(x = plotData[, 2:11],
      pch = 19,  cex = 1,
      col = plotData$ascocarp,
      lower.panel=NULL)


pltTitle <- "Principal Component Analysis"

## decide which PCs to use for plotting
pcToPlot <- c(1, 2)
pcCols <- grep(pattern = "Dim.", x = colnames(plotData), value = T)[pcToPlot]
fillColumn <- "ascocarp"

pointCol <- base::structure(
  .Data = c("green", "red", "blue"),
  names = levels(plotData$ascocarp))

pcaPlot <- ggplot(
  data = plotData,
  mapping = aes(x = !!sym(pcCols[1]), y = !!sym(pcCols[2]), label = gene)) +
  geom_point(mapping = aes(color = !!sym(fillColumn)),
             size = 2, stroke = 2) +
  scale_shape_manual(values = c(1, 15, 17)) +
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
        panel.grid = element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold"))


print(pcaPlot)

png(filename = paste(outPrefix, ".PCA.FPKM.png", sep = ""), width = 5500, height = 5000, res = 550)
print(pcaPlot)
dev.off()


