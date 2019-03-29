library(chipmine)
library(org.AFumigatus293.eg.db)


rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/GO_enrichment/topGO_functions.R")

outDir <- here::here("analysis", "RNAseq_analysis")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

##################################################################################

file_goMap <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_orgDb/geneid2go.AFumigatus_Af293.topGO.map"

orgDb <- org.AFumigatus293.eg.db
##################################################################################


















