library(chipmine)
library(here)
library(BSgenome.Afumigatus.AspGD.Af293)

## this script
## 1) build gene level peak target matrix from curated peak target matrix file

rm(list = ls())

outDir <- here::here("ChIPseq_analysis", "peak_targets")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}


##################################################################################

file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_version_s03-m05-r06_CDS_ORFs.bed"

TF_dataPath <- here::here("data", "TF_data")

sampleList <- c("CREEHA_CONTROL4", "CREEHA_CONTROL5", "CREEHA_10MMAA4", "CREEHA_10MMAA5")
file_targets <- here::here("ChIPseq_analysis", "peak_targets", "peak_targets.common.all.curated.tab")

##################################################################################

## get the sample details
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = sampleList,
                                   dataPath = TF_dataPath,
                                   matrixSource = "normalizedmatrix")

exptDataList <- purrr::transpose(exptData) %>%
  purrr::set_names(nm = purrr::map(., "sampleId"))

tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered", "summitSeq"),
  FUN = function(x){ structure(paste(x, ".", sampleList, sep = ""), names = sampleList) },
  simplify = F, USE.NAMES = T
)

targetsAll <- suppressMessages(readr::read_tsv(file = file_targets, col_names = T))

##################################################################################

## prepare ChIPseq target data
s1 <- "CREEHA_CONTROL4"
s2 <- "CREEHA_10MMAA4"

s1SummitSeq <- get_narrowpeak_summit_seq(npFile = exptDataList[[s1]]$narrowpeakFile,
                                         id = exptDataList[[s1]]$sampleId,
                                         genome = BSgenome.Afumigatus.AspGD.Af293,
                                         length = 200)

s1Targets <- targetsAll %>% 
  dplyr::filter_at(.vars = tfCols$hasPeak[s1], .vars_predicate = all_vars(. == TRUE)) %>%
  dplyr::filter_at(.vars = tfCols$peakPosition[s1], .vars_predicate = all_vars(. == "TSS")) %>%
  dplyr::filter_at(.vars = tfCols$pvalFiltered[s1], .vars_predicate = all_vars(. == TRUE)) %>% 
  dplyr::group_by(gene) %>% 
  dplyr::arrange(desc(!! as.name(tfCols$peakPval[s1])), .by_group = TRUE) %>% 
  dplyr::slice(1L) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(y = s1SummitSeq, by = structure("name", names = tfCols$peakId[s1])) %>% 
  dplyr::select(gene, unname(c(tfCols$hasPeak[s1], tfCols$peakPval[s1], tfCols$peakId[s1], tfCols$summitSeq[s1])))

nrow(s1Targets)
length(unique(s1Targets$gene))


s2SummitSeq <- get_narrowpeak_summit_seq(npFile = exptDataList[[s2]]$narrowpeakFile,
                                         id = exptDataList[[s2]]$sampleId,
                                         genome = BSgenome.Afumigatus.AspGD.Af293,
                                         length = 200)

s2Targets <- targetsAll %>% 
  dplyr::filter_at(.vars = tfCols$hasPeak[s2], .vars_predicate = all_vars(. == TRUE)) %>%
  dplyr::filter_at(.vars = tfCols$peakPosition[s2], .vars_predicate = all_vars(. == "TSS")) %>%
  dplyr::filter_at(.vars = tfCols$pvalFiltered[s2], .vars_predicate = all_vars(. == TRUE)) %>% 
  dplyr::group_by(gene) %>% 
  dplyr::arrange(desc(!! as.name(tfCols$peakPval[s2])), .by_group = TRUE) %>% 
  dplyr::slice(1L) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(y = s2SummitSeq, by = structure("name", names = tfCols$peakId[s2])) %>% 
  dplyr::select(gene, unname(c(tfCols$hasPeak[s2], tfCols$peakPval[s2], tfCols$peakId[s2], tfCols$summitSeq[s2])))

nrow(s2Targets)
length(unique(s2Targets$gene))


targetSet <- data.frame(geneId = union(s1Targets$gene, s2Targets$gene), stringsAsFactors = F) %>% 
  dplyr::left_join(y = s1Targets, by = c("geneId" = "gene")) %>% 
  dplyr::left_join(y = s2Targets, by = c("geneId" = "gene")) %>% 
  tidyr::replace_na(purrr::set_names(list(FALSE, FALSE), unname(tfCols$hasPeak[c(s1, s2)])))


readr::write_tsv(x = targetSet, path = paste(outDir, "/peak_targets.curated.filtered.tab", sep = ""))










