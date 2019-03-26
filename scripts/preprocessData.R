library(chipmine)
library(foreach)
library(doParallel)
library(here)
library(TxDb.Afumigatus.Af293.AspGD.GFF)

rm(list = ls())


##################################################################################

file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_version_s03-m05-r06_CDS_ORFs.bed"
orgDb <- NULL


TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")
other_dataPath <- here::here("data", "other_data")


##################################################################################
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::select(-score) %>% 
  dplyr::mutate(length = end - start)



tfSampleFile <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")

tfSampleList <- readr::read_tsv(file = tfSampleFile, col_names = c("id"),  comment = "#") %>%
  as.data.frame()


tfInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = tfSampleList$id,
                                 dataPath = TF_dataPath,
                                 matrixSource = "normalizedmatrix")


i <- 7

for(i in 1:nrow(tfInfo)){
  
  peakAn <- narrowPeak_annotate(peakFile = tfInfo$narrowpeakFile[i],
                                txdb = TxDb.Afumigatus.Af293.AspGD.GFF,
                                includeFractionCut = 0.7,
                                bindingInGene = FALSE, promoterLength = 500,
                                insideSkewToEndCut = 0.7,
                                output = tfInfo$narrowpeakAnno[i])
  
  
  tfDf <- preProcess_macs2_results(sampleId = tfInfo$sampleId[i],
                                   peakAnnotation = tfInfo$narrowpeakAnno[i],
                                   cdsFile = file_genes,
                                   peakFile = tfInfo$narrowpeakFile[i],
                                   bwFile = tfInfo$bwFile[i],
                                   outFile = tfInfo$tfPeakFile[i],
                                   bindingInGene = FALSE)
  
  
}



##################################################################################




