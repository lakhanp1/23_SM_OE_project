suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))

rm(list = ls())

##################################################################################

file_backbone <- here::here("data", "reference_data", "backbone_enz_type.txt")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

backboneTypes <- suppressMessages(readr::read_tsv(file = file_backbone))

smData <- AnnotationDbi::select(
  x = orgDb, keys = backboneTypes$geneId, columns = c("SM_ID", "SM_CLUSTER"), keytype = "GID"
)

backboneTypes <- dplyr::left_join(x = backboneTypes, y = smData, by = c("geneId" = "GID")) %>% 
  dplyr::arrange(geneId) %>% 
  dplyr::filter(!is.na(SM_ID))

clusterTYpe <- dplyr::group_by(backboneTypes, SM_ID) %>% 
  dplyr::summarise(
    genes = paste(geneId, collapse = ","),
    backbones = paste(sort(unique(backbone_type)), collapse = ",")
  )


clusterBackboneStats <- tidyr::pivot_wider(
  data = backboneTypes,
  id_cols = c(SM_ID, SM_CLUSTER),
  names_from = backbone_type,
  values_from = geneId,
  values_fn = length,
  values_fill = 0
) %>% 
  dplyr::select(
    SM_ID, SM_CLUSTER, PKS, NRPS, terpene, DMATS, PKS_NRPS, everything()
  ) %>% 
  dplyr::arrange(SM_ID)

readr::write_tsv(
  file = here("data", "reference_data", "SM_cluster_backbone.stats.tab"),
  x = clusterBackboneStats
)


clusterBackboneGenes <- tidyr::pivot_wider(
  data = backboneTypes,
  id_cols = c(SM_ID, SM_CLUSTER),
  names_from = backbone_type,
  values_from = geneId,
  values_fn = list(geneId = ~paste(., collapse = ";")),
  values_fill = "0"
) %>% 
  dplyr::select(
    SM_ID, SM_CLUSTER, PKS, NRPS, terpene, DMATS, PKS_NRPS, everything()
  ) %>% 
  dplyr::arrange(SM_ID)


readr::write_tsv(
  file = here("data", "reference_data", "SM_cluster_backbone.genes.tab"),
  x = clusterBackboneGenes
)


