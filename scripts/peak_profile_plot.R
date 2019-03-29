library(chipmine)
library(org.AFumigatus293.eg.db)
library(BSgenome.Afumigatus.AspGD.Af293)
library(DiffBind)


## this script:
## 1) use diffbind regions and generate profile plot around peak


rm(list = ls())

outDir <- here::here("analysis", "ChIPseq_analysis", "diffBind")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

analysisName <- "creE_diffbind"
outPrefix <- paste(outDir, "/", analysisName, sep = "")

##################################################################################

orgDb <- org.AFumigatus293.eg.db

file_sampleInfo <- here::here("analysis", "ChIPseq_analysis", "diffBind", "sampleInfo.txt")
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_targets <- here::here("analysis", "ChIPseq_analysis", "peak_targets", "peak_targets.common.all.curated.tab")
file_diffBindAll <- here::here("analysis", "ChIPseq_analysis", "diffBind", "creE_diffbind.all.tab")
TF_dataPath <- here::here("data", "TF_data")

sampleInfo <- suppressMessages(readr::read_tsv(file = file_exptInfo))

## get the sample details
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = sampleInfo$sampleId,
                                   dataPath = TF_dataPath,
                                   matrixSource = "normalizedmatrix")

exptDataList <- purrr::transpose(exptData) %>%
  purrr::set_names(nm = purrr::map(., "sampleId"))

##################################################################################

# diffDba <- DiffBind::dba.load(file = paste(analysisName, ".dba", sep = ""), dir = outDir, pre = "")

diffbindRes <- suppressMessages(readr::read_tsv(file = file_diffBindAll)) %>% 
  dplyr::select(seqnames, start, end, name, diffBind, peakOccupancy, peakDiff) %>% 
  dplyr::filter(! peakOccupancy %in% c("no_consensus")) %>% 
  dplyr::distinct(name, .keep_all = TRUE)

diffGr <- GenomicRanges::makeGRangesFromDataFrame(diffbindRes, keep.extra.columns = T)

which(table(diffGr$name) > 1)

peakDiffAn <- dplyr::select(diffbindRes, name, peakDiff) %>% 
  dplyr::rename(gene = name, cluster = peakDiff) %>% 
  as.data.frame()

# ## for the first time, generate profile matrices.
# ## next time these matrices can be directly imported
# for(i in 1:nrow(exptData)){
#   mat <- bigwig_profile_matrix(bwFile = exptData$bwFile[i],
#                                regions = diffGr,
#                                signalName = exptData$sampleId[i],
#                                storeLocal = TRUE, localPath = exptData$matFile[i],
#                                extend = c(500, 500), target = "region")
# }



matList <- import_profiles(exptInfo = exptData, geneList = diffGr$name,
                           source = "normalizedmatrix",
                           up = 50, target = 30, down = 50)


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
                      colors = c("white", "red")))
  }
)

ylimList <- sapply(X = exptData$sampleId, FUN = function(x) c(0, 80), simplify = FALSE)

profilePlots <- multi_profile_plots(exptInfo = exptData, genesToPlot = diffGr$name,
                                    profileColors = tfColorList,
                                    clusters = peakDiffAn,
                                    targetType = "region",
                                    matBins = c(50, 30, 50, 10), matSource = "normalizedmatrix",
                                    column_title_gp = gpar(fontsize = 12),
                                    ylimFraction = ylimList)



# pdf(file = paste(outPrefix, "_profiles.pdf", sep = ""), width = 18, height = 13)
png(file = paste(outPrefix, "_profiles.png", sep = ""), width = 4000, height = 2500, res = 250)

ht <- draw(profilePlots$heatmapList,
     main_heatmap = exptData$profileName[1],
     column_title = "creE peaks diffbind comparison",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     heatmap_legend_side = "bottom",
     gap = unit(7, "mm"),
     padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()











