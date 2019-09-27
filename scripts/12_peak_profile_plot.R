library(chipmine)

## this script:
## 1) use diffbind regions and generate profile plot around peak

rm(list = ls())

outDir <- here::here("analysis", "02_ChIPseq_analysis", "03_diffBind")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

analysisName <- "creE_diffbind"
outPrefix <- paste(outDir, "/", analysisName, sep = "")

compare <- c("CREEHA_CONTROL", "CREEHA_10MMAA")

##################################################################################

file_diffbindInfo <- here::here("analysis", "02_ChIPseq_analysis", "03_diffBind", "sampleInfo.txt")
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_diffbindRes <- here::here("analysis", "02_ChIPseq_analysis", "03_diffBind",
                               "creE_diffbind.annotation.filtered.tab")
TF_dataPath <- here::here("data", "TF_data")

sampleInfo <- suppressMessages(readr::read_tsv(file = file_diffbindInfo))

## get the sample details
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   dataPath = TF_dataPath,
                                   profileMatrixSuffix = "normalizedmatrix")

exptDataList <- purrr::transpose(exptData) %>%
  purrr::set_names(nm = purrr::map(., "sampleId"))

grp1 <- compare[1]
grp2 <- compare[2]
grp1Index <- which(exptData$groupId == grp1)
grp2Index <- which(exptData$groupId == grp2)
grp1Samples <- exptData$sampleId[grp1Index]
grp2Samples <- exptData$sampleId[grp2Index]
grp1SpecificOcc = paste(grp1, ":specific", sep = "")
grp2SpecificOcc = paste(grp2, ":specific", sep = "")
grp1EnrichedCategory <- paste(grp1, ":enriched", sep = "")
grp2EnrichedCategory <- paste(grp2, ":enriched", sep = "")

bestGrp1Id <- exptData$sampleId[exptData$groupId == grp1 & exptData$bestRep == 1]
bestGrp2Id <- exptData$sampleId[exptData$groupId == grp2 & exptData$bestRep == 1]

##################################################################################

# diffDba <- DiffBind::dba.load(file = paste(analysisName, ".dba", sep = ""), dir = outDir, pre = "")

diffbindRes <- suppressMessages(readr::read_tsv(file = file_diffbindRes)) %>% 
  dplyr::select(seqnames, start, end, name, Fold, FDR, diffBind, peakOccupancy,
                categoryDiffbind, pvalGood.all) %>% 
  dplyr::filter(pvalGood.all != 0) %>% 
  dplyr::distinct(name, .keep_all = TRUE) %>% 
  dplyr::mutate(
    cluster = dplyr::case_when(
      peakOccupancy == "common" & categoryDiffbind != "common" ~ 
        paste(peakOccupancy, categoryDiffbind, sep = "\n"),
      TRUE ~ categoryDiffbind 
    )
  )

diffGr <- GenomicRanges::makeGRangesFromDataFrame(diffbindRes, keep.extra.columns = T)

##################################################################################
## get the average summit position
peakList <- GenomicRanges::GRangesList(
  lapply(X = exptData$peakFile[c(grp1Index, grp2Index)],
         FUN = rtracklayer::import, format = "narrowPeak")
)

names(peakList) <- exptData$sampleId[c(grp1Index, grp2Index)]

# pgr <- peakList$CREEHA_CONTROL4
## find overlap of each peak GR with diffGr.
## if multiple peaks overlap with a diffGr, select strongest
ovPk <- endoapply(
  X = peakList,
  FUN = function(pgr){
    ovlp <- findOverlaps(query = diffGr, subject = pgr)
    opgr <- pgr[ovlp@to]
    mcols(opgr)$diffGrId <- ovlp@from
    opgr <- as.data.frame(opgr) %>% 
      dplyr::group_by(diffGrId) %>% 
      dplyr::arrange(desc(pValue)) %>% 
      dplyr::slice(1L) %>% 
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    return(opgr)
  })

combinedPeaks <- unlist(ovPk)

summitPos <- GenomicRanges::resize(
  x = GenomicRanges::shift(x = combinedPeaks, shift = combinedPeaks$peak),
  width = 1, fix = "start"
)

avgSummit <- as.data.frame(x = summitPos, row.names = NULL) %>% 
  dplyr::group_by(diffGrId) %>% 
  dplyr::summarise(meanSummit = round(mean(start)),
                   summitSd = round(sd(start))
  )

diffGr$avgSummit[avgSummit$diffGrId] <- avgSummit$meanSummit
diffGr$summitSd[avgSummit$diffGrId] <- avgSummit$summitSd

##################################################################################
## draw summit position variation for each peak region
as.data.frame(diffGr) %>% 
  dplyr::select(name, peakOccupancy, categoryDiffbind, avgSummit, summitSd) %>% 
  dplyr::filter(peakOccupancy == "common") %>% 
  ggplot(mapping = aes(x = summitSd)) +
  geom_histogram(binwidth = 5) +
  scale_x_continuous(breaks = seq(0, 500, by = 50)) +
  labs(title = "peak summit standard deviation for common peaks") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold")
  )

summitPos$sampleId <- names(summitPos)

summitSummaryDf <- as.data.frame(summitPos, row.names = NULL) %>% 
  dplyr::select(name, start, diffGrId, sampleId) %>% 
  tidyr::pivot_wider(id_cols = diffGrId, names_from = sampleId, values_from = start) %>% 
  dplyr::mutate(AA_CTRL4 = abs(CREEHA_10MMAA4 - CREEHA_CONTROL4),
                AA_CTRL5 = abs(CREEHA_10MMAA5 - CREEHA_CONTROL5)) %>% 
  dplyr::mutate(AA_CTRL_avg = (AA_CTRL4 + AA_CTRL5) / 2)

summitSummaryDf$peakOccupancy <- diffGr$peakOccupancy[summitSummaryDf$diffGrId]
summitSummaryDf$categoryDiffbind <- diffGr$categoryDiffbind[summitSummaryDf$diffGrId]
summitSummaryDf$name <- diffGr$name[summitSummaryDf$diffGrId]

summitSummaryDf <- dplyr::select(
  summitSummaryDf, diffGrId, name, peakOccupancy, categoryDiffbind,
  starts_with("AA_CTRL"))

dplyr::filter_at(.tbl = summitSummaryDf, 
                 .vars = vars(starts_with("AA_CTRL")),
                 .vars_predicate = any_vars(. > 50 & . < 100)) %>% 
  dplyr::arrange(desc(AA_CTRL_avg)) %>% 
  dplyr::filter(peakOccupancy == "common")


dplyr::filter(summitSummaryDf, peakOccupancy == "common") %>% 
  ggplot() +
  geom_histogram(mapping = aes(x = AA_CTRL4), fill = "red", binwidth = 5, alpha = 0.4) +
  geom_histogram(mapping = aes(x = AA_CTRL5), fill = "blue", binwidth = 5, alpha = 0.4) +
  scale_x_continuous(breaks = seq(0, 1000, by = 100)) +
  labs(title = "peak summit difference: common peaks",
       x = "summit position difference (bp)") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 14, face = "bold")
  )


##################################################################################

# peakCenterGr <- GenomicRanges::resize(x = diffGr, fix = "center", width = 1, use.names = T)
peakCenterGr <- GRanges(seqnames = seqnames(diffGr),
                        ranges = IRanges(start = diffGr$avgSummit, width = 1))

mcols(peakCenterGr) <- mcols(diffGr)

which(table(diffGr$name) > 1)

peakDiffAn <- dplyr::select(diffbindRes, name, Fold, FDR, diffBind, peakOccupancy,
                            categoryDiffbind, pvalGood.all, cluster) %>% 
  dplyr::mutate(
    pvalGroup = if_else(pvalGood.all == 0, "weak", "strong"),
    rankMetric = (-log10(FDR) * sign(Fold))) %>%
  dplyr::rename(geneId = name) %>%
  as.data.frame()

peakDiffAn$diffBind <- factor(x = peakDiffAn$diffBind, levels = c("down", "noDiff", "up"))
peakDiffAn$peakOccupancy <- factor(
  x = peakDiffAn$peakOccupancy,
  levels = c(grp1SpecificOcc, "common", grp2SpecificOcc))
peakDiffAn$cluster <- factor(
  x = peakDiffAn$cluster,
  levels = c(grp1SpecificOcc, paste("common", grp1EnrichedCategory, sep = "\n"),
             "common",
             paste("common", grp2EnrichedCategory, sep = "\n"), grp2SpecificOcc))

# ## for the first time, generate profile matrices.
# ## next time these matrices can be directly imported
# for(i in 1:nrow(exptData)){
#   mat <- bigwig_profile_matrix(bwFile = exptData$bwFile[i],
#                                regions = peakCenterGr,
#                                signalName = exptData$sampleId[i],
#                                storeLocal = TRUE, localPath = exptData$matFile[i],
#                                extend = c(1000, 1000), targetName = "peak center")
# }


matList <- import_profiles(exptInfo = exptData, geneList = diffGr$name,
                           source = "normalizedmatrix",
                           targetType = "point", targetName = "peak center",
                           up = 100, target = 0, down = 100)


## tf colors
tfMeanProfile <- NULL
if(length(exptData$sampleId) == 1){
  tfMeanProfile <- matList[[exptData$sampleId[1]]]
} else{
  tfMeanProfile <- getSignalsFromList(lt = matList[exptData$sampleId])
}

quantile(tfMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# tfMeanColor <- colorRamp2(quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T), c("white", "red"))
tfColorList <- sapply(
  X = exptData$sampleId,
  FUN = function(x){
    return(colorRamp2(breaks = quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T),
                      colors = unlist(strsplit(x = exptDataList[[x]]$color, split = ",")))
    )
  }
)

# ylimList <- sapply(X = exptData$sampleId, FUN = function(x) c(0, 80), simplify = FALSE)
ylimList <- list(
  CREEHA_CONTROL4 = c(0, 75), CREEHA_CONTROL5 = c(0, 55), WT_CONTROL5 = c(0, 70),
  CREEHA_10MMAA4 = c(0, 115), CREEHA_10MMAA5 = c(0, 115), WT_10MMAA5 = c(0, 70)
)

profilePlots <- multi_profile_plots(
  exptInfo = exptData, genesToPlot = diffGr$name,
  profileColors = tfColorList,
  clusters = dplyr::select(peakDiffAn, geneId, cluster),
  showAnnotation = FALSE,
  clustOrd = levels(peakDiffAn$cluster),
  targetType = "point",
  targetName = "summit",
  matBins = c(100, 0, 100, 10), matSource = "normalizedmatrix",
  column_title_gp = gpar(fontsize = 12),
  ylimFraction = ylimList
)



# pdf(file = paste(outPrefix, "_profiles.pdf", sep = ""), width = 18, height = 13)
png(file = paste(outPrefix, "_profiles.png", sep = ""), width = 3500, height = 2500, res = 250)

ht <- draw(
  profilePlots$heatmapList,
  main_heatmap = exptData$profileName[1],
  # split = dplyr::select(peakDiffAn, pvalGroup, diffBind),
  column_title = "creE peaks diffbind comparison",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_sub_title_side = "left",
  heatmap_legend_side = "bottom",
  gap = unit(7, "mm"),
  padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()




htList2 <- profilePlots$profileHeatmaps$CREEHA_CONTROL4$heatmap +
  profilePlots$profileHeatmaps$CREEHA_10MMAA4$heatmap +
  profilePlots$profileHeatmaps$WT_CONTROL5$heatmap +
  profilePlots$profileHeatmaps$WT_10MMAA5$heatmap

# pdf(file = paste(outPrefix, "_profiles.pdf", sep = ""), width = 18, height = 13)
png(file = paste(outPrefix, ".profiles.ungrouped.png", sep = ""), width = 3000, height = 3000, res = 250)

ht <- draw(
  htList2,
  # main_heatmap = exptData$profileName[1],
  split = rep(1, nrow(tfMeanProfile)),
  column_title = "creE peaks diffbind comparison",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_sub_title_side = "left",
  heatmap_legend_side = "bottom",
  gap = unit(7, "mm"),
  padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()

##################################################################################


## TF1 scalled matrix
tf1ScalledMat <- scale(x = matList[[bestGrp1Id]], center = TRUE, scale = TRUE)

quantile(tf1ScalledMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
tf1ScalledColor <- colorRamp2(
  breaks = quantile(tf1ScalledMat, c(0.50, 0.99), na.rm = T),
  colors = unlist(stringr::str_split(exptDataList[[bestGrp1Id]]$color, pattern = ","))
)

## TF2 scalled matrix
tf2ScalledMat <- scale(x = matList[[bestGrp2Id]], center = TRUE, scale = TRUE)
quantile(tf2ScalledMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
tf2ScalledColor <- colorRamp2(
  breaks = quantile(tf2ScalledMat, c(0.50, 0.99), na.rm = T),
  colors = unlist(stringr::str_split(exptDataList[[bestGrp2Id]]$color, pattern = ","))
)

## Difference between TF2 and TF1 scalled matrix
scalledTfDiffMat <- tf2ScalledMat - tf1ScalledMat
attr(scalledTfDiffMat, "signal_name") <- "fold_change"
attr(scalledTfDiffMat, "target_name") <- "summit"
plot(density(scalledTfDiffMat))
quantile(scalledTfDiffMat, c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)


scalledTfDiffColor <- colorRamp2(breaks = c(-3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3),
                                 colors = RColorBrewer::brewer.pal(n = 11, name = "PRGn"))



## Scalled TF diff profile
scalledTfDiffProf <- profile_heatmap(
  profileMat = scalledTfDiffMat,
  showAnnotation = FALSE,
  signalName = "fold_change",
  profileColor = scalledTfDiffColor,
  column_title_gp = gpar(fontsize = 12),
  ylimFraction = c(-2.3, 1.7)
)

rowOrd <- order(peakDiffAn$rankMetric, decreasing = TRUE)

htListDiff <- scalledTfDiffProf$heatmap +
  profilePlots$profileHeatmaps[[bestGrp2Id]]$heatmap +
  profilePlots$profileHeatmaps[[bestGrp1Id]]$heatmap


pdf(file = paste(outPrefix, ".diff_profiles.pdf", sep = ""), width = 10, height = 12)
# png(file = paste(outPrefix, ".profiles.ungrouped.png", sep = ""), width = 3000, height = 3000, res = 250)

ht <- draw(
  htListDiff,
  main_heatmap = "fold_change",
  # row_order = rowOrd,
  column_title = paste("rglT peaks ", bestGrp2Id, "/", bestGrp1Id, " signal ratio", sep = ""),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_sub_title_side = "left",
  heatmap_legend_side = "bottom",
  gap = unit(7, "mm"),
  padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()




