suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))

rm(list = ls())

##################################################################################

file_exptInfo <- here::here("data", "referenceData/sample_info.txt")

file_genes <- here::here("data", "referenceData/AN_genes_for_polII.bed")

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")
other_dataPath <- here::here("data", "other_data")

file_tf_macs2 <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")
file_tf <- paste(TF_dataPath, "/", "sample_tf.list", sep = "")


set.seed(20)
clusters <- 7
tfYlim <- 0.996              ##0.999
geneFilter <- c("AN5245", "AN3245")

cl <- makeCluster(3) #not to overload your computer
registerDoParallel(cl)

##################################################################################

geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "name", "score", "strand")) %>% 
  dplyr::mutate(length = end - start) %>% 
  dplyr::filter(! name %in% geneFilter)


tfSampleList <- readr::read_tsv(file = file_tf_macs2, col_names = c("id"),  comment = "#")

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSampleList$id,
  dataPath = TF_dataPath,
  matrixSource = "normalizedmatrix")

i <- 1

foreach(i = 1:nrow(tfInfo),
        .packages = c("chipmine")) %dopar% {
          
          ## 2kb - 2kb - 1kb matrix
          bwMat <- chipmine::bigwig_profile_matrix(bwFile = tfInfo$bwFile[i],
                                                   regions = file_genes,
                                                   signalName = tfInfo$sampleId[i],
                                                   extend = c(2000, 1000), w = 10,
                                                   storeLocal = TRUE,
                                                   localPath = tfInfo$matFile[i],
                                                   target_ratio = 0.4)
          
          mat1 <- chipmine::import_profile_from_file(
            file = tfInfo$matFile[i],
            source = "normalizedmatrix",
            signalName = tfInfo$sampleId[i],
            selectGenes = geneSet$name)
          
          ## check the distribution in data
          quantile(mat1, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
          # col_fun <- colorRamp2(quantile(mat1, c(0.50, 0.995), na.rm = T), c("white", "red"))
          
          km <- chipmine::profile_matrix_kmeans(
            mat = mat1,
            km = clusters,
            clustFile = tfInfo$clusterFile[i],
            name = tfInfo$sampleId[i])
          
          
          cat(as.character(Sys.time()), "Done...", tfInfo$sampleId[i], "\n")
          
        }


stopCluster(cl)








