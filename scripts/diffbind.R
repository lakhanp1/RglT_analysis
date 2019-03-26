library(chipmine)
library(org.AFumigatus293.eg.db)
library(BSgenome.Afumigatus.AspGD.Af293)
library(DiffBind)

## this script:
## 1) run DiffBind to perform differential binding analysis
## 2) add peak data for individual sample
## 3) prepare combined report
## 4) add peak target gene information for each best peaksets

rm(list = ls())

outDir <- here::here("ChIPseq_analysis", "diffBind")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

compare <- c("CREEHA_CONTROL", "CREEHA_10MMAA")

analysisName <- "creE_diffbind"
outPrefix <- paste(outDir, "/", analysisName, sep = "")

fdr_cut <- 0.05
lfc_cut <- 1
up_cut <- lfc_cut
down_cut <- lfc_cut * -1
##################################################################################
orgDb <- org.AFumigatus293.eg.db

file_sampleInfo <- here::here("ChIPseq_analysis", "diffBind", "sampleInfo.txt")
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_targets <- here::here("ChIPseq_analysis", "peak_targets", "peak_targets.common.all.curated.tab")

TF_dataPath <- here::here("data", "TF_data")

##################################################################################
sampleInfo <- suppressMessages(readr::read_tsv(file = file_sampleInfo))

## get the sample details
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = sampleInfo$SampleID,
                                   dataPath = TF_dataPath,
                                   matrixSource = "normalizedmatrix")

exptDataList <- purrr::transpose(exptData) %>%
  purrr::set_names(nm = purrr::map(., "sampleId"))

tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered", "summitSeq", "overlap", "targetGene"),
  FUN = function(x){ structure(paste(x, ".", exptData$sampleId, sep = ""), names = exptData$sampleId) },
  simplify = F, USE.NAMES = T
)

grp1 <- compare[1]
grp2 <- compare[2]
grp1Index <- which(exptData$groupId == grp1)
grp2Index <- which(exptData$groupId == grp2)
grp1Samples <- exptData$sampleId[grp1Index]
grp2Samples <- exptData$sampleId[grp2Index]

groupCols <- sapply(
  X = c("peakCall"),
  FUN = function(x){ structure(paste(x, ".", compare, sep = ""), names = compare) },
  simplify = F, USE.NAMES = T
)

##################################################################################
## diffBind analysis

dataDba <- DiffBind::dba(sampleSheet = sampleInfo)

plot(dataDba)

countsDba <- DiffBind::dba.count(DBA = dataDba, score = "DBA_SCORE_TMM_MINUS_FULL", summits = 150)

plot(countsDba)

countsDba <- DiffBind::dba.contrast(DBA = countsDba, categories = DBA_CONDITION, minMembers = 2)

diffDba <- DiffBind::dba.analyze(DBA = countsDba, method = DBA_DESEQ2, bReduceObjects = FALSE)
DiffBind::dba.save(DBA = diffDba, file = paste(analysisName, ".dba", sep = ""), dir = outDir, pre = "")

plot(diffDba, contrast=1)

diffDf <- DiffBind::dba.report(DBA = diffDba, bFlip = TRUE, th = 1,
                               bCalledDetail = T, DataType = DBA_DATA_FRAME)

diffDf <- dplyr::mutate(diffDf, region = paste(Chr, ":", Start, "-", End, sep = ""))

# dba.plotMA(diffDba)
# dba.plotVolcano(diffDba)
# dba.plotBox(diffDba)
# dba.plotHeatmap(diffDba, contrast=1, correlations=FALSE)
# 
# dba.overlap(DBA = diffDba, mask = diffDba$masks$All, mode = DBA_OLAP_PEAKS)
# dba.peakset(DBA = diffDba, consensus = )

##################################################################################
## create report table

## add individual peak information columns
diffGr <- DiffBind::dba.report(DBA = diffDba, bFlip = TRUE, th = 1)

## add overlapping peaks from each sample
diffRes <- combinatorial_binding_matrix(sampleInfo = exptData, peakRegions = diffGr)


## count of samples which showed macs2 peak under condition1
diffRes[[unname(groupCols$peakCall[grp1])]] <- purrr::pmap_int(
  .l = dplyr::select(diffRes, unname(tfCols$overlap[grp1Samples])),
  .f = sum, na.rm = TRUE
)

## count of samples which showed macs2 peak under condition2
diffRes[[unname(groupCols$peakCall[grp2])]] <- purrr::pmap_int(
  .l = dplyr::select(diffRes, unname(tfCols$overlap[grp2Samples])),
  .f = sum, na.rm = TRUE
)


dplyr::group_by_at(diffRes, .vars = vars(starts_with("peakCall"))) %>% 
  dplyr::summarise(n = n())

## add diffBind status and peak occupancy status
diffData <- diffRes %>% 
  dplyr::mutate(
    diffBind = dplyr::case_when(
      Fold >= up_cut & FDR <= fdr_cut ~ "up",
      Fold <= down_cut & FDR <= fdr_cut ~ "down",
      TRUE ~ "noDiff"
    ),
    peakOccupancy = dplyr::case_when(
      !!as.name(groupCols$peakCall[grp1]) >= 2 & !!as.name(groupCols$peakCall[grp2]) >= 2 ~ "common",
      !!as.name(groupCols$peakCall[grp1]) >= 2 & !!as.name(groupCols$peakCall[grp2]) < 2 ~ 
        paste("specific:", grp1, sep = ""),
      !!as.name(groupCols$peakCall[grp1]) < 2 & !!as.name(groupCols$peakCall[grp2]) >= 2 ~ 
        paste("specific:", grp2, sep = ""),
      TRUE ~ "no_consensus"
    )
  )


dplyr::group_by(diffData, peakOccupancy, diffBind) %>% 
  dplyr::summarise(n = n())

## final diffbind data
diffData <- dplyr::mutate(diffData, peakDiff = paste(peakOccupancy, diffBind, sep = ":")) %>% 
  dplyr::select(seqnames, start, end, name, starts_with("Conc"), Fold, p.value, FDR,
                diffBind, peakOccupancy, peakDiff, starts_with("peakCall."), everything())


readr::write_tsv(x = diffData, path = paste(outPrefix, "_all.tab", sep = ""))

##################################################################################

bestGrp1Id <- exptData$sampleId[exptData$groupId == grp1 & exptData$bestRep == 1]
bestGrp2Id <- exptData$sampleId[exptData$groupId == grp2 & exptData$bestRep == 1]

targetsAll <- suppressMessages(readr::read_tsv(file = file_targets, col_names = T))

## extract targets specific for each sample
s1Targets <- dplyr::select(targetsAll, gene, contains(bestGrp1Id)) %>% 
  dplyr::filter_at(.vars = tfCols$hasPeak[bestGrp1Id], .vars_predicate = all_vars(. == TRUE)) %>% 
  dplyr::select(gene, starts_with("peakId"), starts_with("peakPosition"),
                starts_with("peakType"), starts_with("peakDist"), starts_with("pvalFiltered")) %>% 
  dplyr::distinct() %>% 
  dplyr::rename(!! as.name(tfCols$targetGene[bestGrp1Id]) := gene)

s2Targets <- dplyr::select(targetsAll, gene, contains(bestGrp2Id)) %>% 
  dplyr::filter_at(.vars = tfCols$hasPeak[bestGrp2Id], .vars_predicate = all_vars(. == TRUE)) %>% 
  dplyr::select(gene, starts_with("peakId"), starts_with("peakPosition"),
                starts_with("peakType"), starts_with("peakDist"), starts_with("pvalFiltered")) %>% 
  dplyr::distinct() %>% 
  dplyr::rename(!! as.name(tfCols$targetGene[bestGrp2Id]) := gene)


diffAnn <- diffData %>% 
  dplyr::arrange(seqnames, start) %>% 
  dplyr::left_join(y = s1Targets, by = unname(tfCols$peakId[bestGrp1Id])) %>% 
  dplyr::left_join(y = s2Targets, by = unname(tfCols$peakId[bestGrp2Id])) %>% 
  dplyr::select(seqnames, start, end, name, starts_with("Conc_"), Fold, p.value, FDR, diffBind, peakOccupancy,
                peakDiff, starts_with("peakCall."), contains(bestGrp1Id), contains(bestGrp2Id)) %>% 
  dplyr::filter((!peakDiff %in% c("specific:CREEHA_CONTROL:noDiff", "specific:CREEHA_10MMAA:noDiff")) &
                  peakOccupancy != "no_consensus") %>% 
  dplyr::mutate(
    targetMatch = if_else(!! as.name(tfCols$targetGene[bestGrp1Id]) == !! as.name(tfCols$targetGene[bestGrp2Id]),
                          TRUE, FALSE)
  ) %>% 
  dplyr::mutate(
    geneId = dplyr::case_when(
      targetMatch ~ !! as.name(tfCols$targetGene[bestGrp1Id]),
      is.na(!! as.name(tfCols$targetGene[bestGrp1Id])) & is.na(!! as.name(tfCols$targetGene[bestGrp2Id])) ~ "NA",
      is.na(!! as.name(tfCols$targetGene[bestGrp1Id])) ~ !! as.name(tfCols$targetGene[bestGrp2Id]),
      is.na(!! as.name(tfCols$targetGene[bestGrp2Id])) ~ !! as.name(tfCols$targetGene[bestGrp1Id]),
      TRUE ~ "NA"
  )) %>% 
  dplyr::select(seqnames, start, end, name, geneId, everything())


diffAnn$pvalFilteredN <- purrr::pmap_int(
  .l = dplyr::select(diffAnn, starts_with("pvalFiltered")),
  .f = sum, na.rm = TRUE
)

readr::write_tsv(x = diffAnn, path = paste(outPrefix, ".annotation.tab", sep = ""))




