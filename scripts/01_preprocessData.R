library(chipmine)
library(foreach)
library(doParallel)
library(here)
library(TxDb.Afumigatus.Af293.AspGD.GFF)
library(org.AFumigatus.Af293.eg.db)

rm(list = ls())


##################################################################################

file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_version_s03-m05-r06_CDS_ORFs.bed"
orgDb <- org.AFumigatus.Af293.eg.db
txDb <- TxDb.Afumigatus.Af293.AspGD.GFF

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")
other_dataPath <- here::here("data", "other_data")


##################################################################################

## prepare txIds excluding rRNA and tRNA transcripts
geneInfo <- AnnotationDbi::select(x = orgDb,
                                  keys = keys(orgDb, keytype = "GID"),
                                  columns = c("GENE_NAME", "DESCRIPTION", "TYPE"),
                                  keytype = "GID") %>% 
  dplyr::rename(geneId = GID)

geneInfo %>% dplyr::filter(!grepl(pattern = "ORF\\|", x = TYPE, perl = TRUE)) %>% 
  dplyr::select(geneId, GENE_NAME, TYPE) %>%
  dplyr::count(TYPE)

geneInfo <- dplyr::filter(geneInfo, !grepl(pattern = "(rRNA|tRNA)\\|", x = TYPE, perl = TRUE))

txInfo <- AnnotationDbi::select(x = txDb, keys = geneInfo$geneId,
                                columns = "TXID", keytype = "GENEID")

tfSampleFile <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")

tfSampleList <- readr::read_tsv(file = tfSampleFile, col_names = c("id"),  comment = "#") %>%
  as.data.frame()


tfInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = tfSampleList$id,
                                 dataPath = TF_dataPath,
                                 profileMatrixSuffix = "normalizedmatrix")


i <- 1

for(i in 1:nrow(tfInfo)){

peakAn <- narrowPeak_annotate(peakFile = tfInfo$peakFile[i],
                              txdb = txDb, promoterLength = 500,
                              txIds = txInfo$TXID,
                              includeFractionCut = 0.7,
                              bindingInGene = FALSE,
                              insideSkewToEndCut = 0.7,
                              output = tfInfo$peakAnno[i],
                              removePseudo = TRUE)


# tfDf <- gene_level_peak_annotation(
#   sampleId = tfInfo$sampleId[i],
#   peakAnnotation = tfInfo$peakAnno[i],
#   preference = c("nearStart", "peakInFeature", "featureInPeak", "upstreamTss", "nearEnd", "intergenic"),
#   genesDf = geneInfo,
#   peakFile = tfInfo$peakFile[i],
#   bwFile = tfInfo$bwFile[i],
#   outFile = tfInfo$peakTargetFile[i])


}



##################################################################################




