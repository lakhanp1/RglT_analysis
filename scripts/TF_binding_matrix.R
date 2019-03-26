library(chipmine)
library(BSgenome.Afumigatus.AspGD.Af293)
library(here)


## This script generate a matrix of all TF 

rm(list = ls())


outDir <- here::here("ChIPseq_analysis", "peak_targets")


if(!dir.exists(outDir)){
  dir.create(path = outDir)
}


file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_version_s03-m05-r06_CDS_ORFs.bed"


TF_dataPath <- here::here("data", "TF_data")

# sampleList <- c("CREEHA_CONTROL2", "CREEHA_CONTROL3", "CREEHA_CONTROL4", "CREEHA_CONTROL5",
#                 "CREEHA_10MMAA2", "CREEHA_10MMAA3", "CREEHA_10MMAA4", "CREEHA_10MMAA5")
# outPrefix <- paste(outDir, "/all_samples", sep = "")

sampleList <- c("CREEHA_CONTROL4", "CREEHA_CONTROL5", "CREEHA_10MMAA4", "CREEHA_10MMAA5")
outPrefix <- paste(outDir, "/good_samples", sep = "")

tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakCoverage",
        "pvalFiltered"),
  FUN = function(x){ structure(paste(x, ".", sampleList, sep = ""), names = sampleList) },
  simplify = F, USE.NAMES = T
)

##################################################################################

genesDf <- data.table::fread(file = file_genes, sep = "\t", header = F, select = c(4), col.names = c("gene"), stringsAsFactors = F)

## get the sample details
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = sampleList,
                                   dataPath = TF_dataPath,
                                   matrixSource = "normalizedmatrix")

exptDataList <- purrr::transpose(exptData) %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

##################################################################################
## generate TF binding matrix based on the target genes

tfData <- get_TF_binding_data(genesDf = genesDf, exptInfo = exptData) %>%
  dplyr::select(gene, starts_with("hasPeak"), starts_with("peakRegion"), starts_with("peakType"),
                starts_with("peakId"), starts_with("peakPval"), starts_with("peakEnrichment"), starts_with("peakDist")) %>%
  dplyr::filter_at(.vars = vars(starts_with("peakType")),
                   .vars_predicate = any_vars(!is.na(.)))

# "hasPeak", "peakPosition", "peakType", "peakId", "peakEnrichment", "peakPval", "peakQval",
# "peakSummit", "peakDist", "summitDist", "bidirectional", "featureCovFrac", "relativeSummitPos",
# "peakRegion", "peakCoverage"

readr::write_tsv(x = tfData, path = paste(outPrefix, "_peak_target_matrix.tab", sep = ""))


##################################################################################
## merge the peaks and generate the combination matrix for presence of peak

## generate the combinatorial binding matrix
mat <- combinatorial_binding_matrix(sampleInfo = exptData,
                                    genome = BSgenome.Afumigatus.AspGD.Af293,
                                    summitSeqLen = 200)

sortCols <- grep(pattern = "^peakPval", x = names(mat), value = T, perl = T)


newMat <- dplyr::group_by_at(.tbl = mat, .vars = vars(starts_with("overlap"))) %>% 
  # dplyr::arrange(desc(`peakPval.CREEHA_CONTROL2`), desc(`peakPval.CREEHA_10MMAA3`), .by_group = TRUE) %>% 
  dplyr::ungroup()


readr::write_tsv(x = newMat, path = paste(outPrefix, "_peak_combinatorial_mat.tab", sep = ""))


# "peakPval.CREEHA_CONTROL2", "peakPval.CREEHA_10MMAA3" 
ggplot2::ggplot(data = mat, mapping = aes(x = `peakPval.CREEHA_10MMAA5`)) +
  geom_histogram(binwidth = 2) +
  coord_cartesian(xlim = c(0, 80))


##################################################################################
## confident peak target matrix based on replicates

## pval cutoff: set to 1 to include all the peaks. else cutoff values will be used from sampleInfo 
# exptData$pval_cutoff <- 1

## get the peaks from best sample among replicates for each condition
s1Id <- exptData[which(exptData$bestRep == 1 & exptData$groupId == "CREEHA_CONTROL"), ]
s1Df <- best_replicate_peakset(sampleInfo = dplyr::filter(exptData, groupId == "CREEHA_CONTROL"),
                               cdsFile = file_genes) %>% 
  dplyr::mutate(
    !! unname(tfCols$pvalFiltered[s1Id$sampleId]) := if_else(
      condition = !! as.name(tfCols$peakPval[s1Id$sampleId]) >= s1Id$pval_cutoff,
      true = TRUE, false = FALSE
    )
  )


masterSet <- dplyr::left_join(x = genesDf, y = s1Df, by = c("gene" = "gene"))

s2Id <- exptData[which(exptData$bestRep == 1 & exptData$groupId == "CREEHA_10MMAA"), ]
s2Df <- best_replicate_peakset(sampleInfo = dplyr::filter(exptData, groupId == "CREEHA_10MMAA"),
                               cdsFile = file_genes) %>% 
  dplyr::mutate(
    !! unname(tfCols$pvalFiltered[s2Id$sampleId]) := if_else(
      condition = !! as.name(tfCols$peakPval[s2Id$sampleId]) >= s2Id$pval_cutoff,
      true = TRUE, false = FALSE
    )
  )

masterSet <- dplyr::left_join(x = masterSet, y = s2Df, by = c("gene" = "gene"))

targets <- dplyr::filter_at(.tbl = masterSet, .vars = vars(starts_with("hasPeak.")),
                            .vars_predicate = any_vars(. == TRUE)) %>% 
  tidyr::replace_na(
    replace = set_names(list(FALSE, FALSE),
                        nm = unname(tfCols$hasPeak[exptData$sampleId[which(exptData$bestRep == 1)]]))) %>% 
  # dplyr::filter_at(.vars = vars(starts_with("pvalFiltered.")), .vars_predicate = any_vars(. == TRUE)) %>% 
  dplyr::arrange(gene)

dplyr::group_by_at(.tbl = targets, .vars = vars(starts_with("hasPeak."))) %>% 
  dplyr::summarise(n = n())

readr::write_tsv(x = targets, path = paste(outDir, "/peak_targets.common.all.tab", sep = ""))

##################################################################################


npGr <- rtracklayer::import(con = exptDataList$CREEHA_10MMAA5$narrowpeakFile, format = "narrowPeak")


quantile(npGr$pValue, probs = c(0, 0.05, seq(0.1, 0.9, 0.1), 0.95, 1))
table(cut(x = npGr$pValue,
          breaks = quantile(npGr$pValue, probs = c(0, 0.05, seq(0.1, 0.9, 0.1), 0.95, 1)),
          labels = paste(c(0, 0.05, seq(0.1, 0.9, 0.1), 0.95) * 100, "-", c(0.05, seq(0.1, 0.9, 0.1), 0.95, 1) * 100, "%", sep = ""),
          right=FALSE))


purrr::map_dfr(
  .x = exptDataList,
  .f = function(f = .x){
    npGr <- rtracklayer::import(con = f$narrowpeakFile, format = "narrowPeak")
    quantile(npGr$pValue, probs = c(0, 0.05, seq(0.1, 0.9, 0.1), 0.95, 1))
    # table(cut(x = npGr$pValue,
    #           breaks = quantile(npGr$pValue, probs = c(0, 0.05, seq(0.1, 0.9, 0.1), 0.95, 1)),
    #           labels = paste(c(0, 0.05, seq(0.1, 0.9, 0.1), 0.95) * 100, "-",
    #                          c(0.05, seq(0.1, 0.9, 0.1), 0.95, 1) * 100, "%", sep = ""),
    #           right=FALSE))
  }) %>% 
  dplyr::mutate(quantile = paste(c(0, 0.05, seq(0.1, 0.9, 0.1), 0.95, 1) * 100, "%", sep = ""))


purrr::map_dfr(
  .x = exptDataList,
  .f = function(f = .x){
    npGr <- rtracklayer::import(con = f$narrowpeakFile, format = "narrowPeak")
    quantile(npGr$pValue, probs = c(0, 0.05, seq(0.1, 0.9, 0.1), 0.95, 1))
  }) %>% 
  dplyr::mutate(quantile = paste(c(0, 0.05, seq(0.1, 0.9, 0.1), 0.95, 1) * 100, "%", sep = ""))











