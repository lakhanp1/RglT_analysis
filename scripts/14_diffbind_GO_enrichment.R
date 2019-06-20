library(chipmine)
library(org.AFumigatus293.eg.db)
library(cowplot)

rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/GO_enrichment/topGO_functions.R")

outDir <- here::here("analysis", "ChIPseq_analysis", "diffbind", "topGO")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

outPrefix <- paste(outDir, "/topGO_goodPeaks.", sep = "")

##################################################################################
diffbindCompare <- c("CREEHA_CONTROL", "CREEHA_10MMAA")

file_goMap <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_orgDb/geneid2go.AFumigatus_Af293.topGO.map"
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")

TF_dataPath <- here::here("data", "TF_data")

sampleList <- c("CREEHA_CONTROL4", "CREEHA_CONTROL5", "CREEHA_10MMAA4", "CREEHA_10MMAA5")
file_diffbindTargets <- here::here("analysis", "ChIPseq_analysis",
                                   "peak_targets", "diffbind_allPeak_targets.tab")
orgDb <- org.AFumigatus293.eg.db

##################################################################################

grp1 <- diffbindCompare[1]
grp2 <- diffbindCompare[2]
grp1Enrich = paste(grp1, ":enriched", sep = "")
grp2Enrich = paste(grp2, ":enriched", sep = "")
grp1Specific = paste(grp1, ":specific", sep = "")
grp2Specific = paste(grp2, ":specific", sep = "")


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
        "peakCoverage", "pvalFiltered"),
  FUN = function(x){ structure(paste(x, ".", sampleList, sep = ""), names = sampleList) },
  simplify = F, USE.NAMES = T
)

diffbindTargets <- suppressMessages(readr::read_tsv(file = file_diffbindTargets, col_names = T))
diffbindTargets <- dplyr::filter(diffbindTargets, pvalFilteredN > 0)

dplyr::group_by(diffbindTargets, categoryDiffbind) %>% 
  dplyr::summarise(n = n())

##################################################################################

## genes specific to CREEHA_CONTROL sample
nodeSize <- 5
outPrefixTf1 <- paste(outPrefix, grp1, "_enriched", sep = "")

tf1StrongTargets <- dplyr::filter(.data = diffbindTargets,
                                  categoryDiffbind %in% c(grp1Enrich, grp1Specific))

nrow(tf1StrongTargets)

tf1Go <- topGO_enrichment(genes = unique(tf1StrongTargets$geneId),
                          goMapFile = file_goMap, goNodeSize = nodeSize)

tf1Go$condition <- grp1Enrich

readr::write_tsv(x = tf1Go, path = paste(outPrefixTf1, ".tab", sep = ""))

plotTitle <- paste(grp1, ": increased binding signal ( n =", length(unique(tf1StrongTargets$geneId)), ")")

tf1GoScatter <- topGO_scatterPlot(df = tf1Go, title = plotTitle)

##################################################################################
## genes specific to CREEHA_10MMAA3 sample

outPrefixTf2 <- paste(outPrefix, grp2, "_enriched", sep = "")

tf2StrongTargets <- dplyr::filter(.data = diffbindTargets,
                                  categoryDiffbind %in% c(grp2Enrich, grp2Specific))
nrow(tf2StrongTargets)

tf2Go <- topGO_enrichment(genes = unique(tf2StrongTargets$geneId),
                          goMapFile = file_goMap, goNodeSize = nodeSize)

tf2Go$condition <- grp2Enrich

readr::write_tsv(x = tf2Go, path = paste(outPrefixTf2, ".tab", sep = ""))

plotTitle <- paste(grp2, ": increased binding signal ( n =", length(unique(tf2StrongTargets$geneId)), ")")
tf2GoScatter <- topGO_scatterPlot(df = tf2Go, title = plotTitle)

##################################################################################
## common targets for CREEHA_CONTROL2 and CREEHA_10MMAA3 samples

outPrefixCommon <- paste(outPrefix, "common_noDiff", sep = "")

nodiffTargets <- dplyr::filter(.data = diffbindTargets, categoryDiffbind == "common")

nrow(nodiffTargets)

commonGo <- topGO_enrichment(genes = unique(nodiffTargets$geneId),
                             goMapFile = file_goMap, goNodeSize = nodeSize)

commonGo$condition <- "common"

readr::write_tsv(x = commonGo, path = paste(outPrefixCommon, ".tab", sep = ""))

plotTitle <- paste("no change in binding signal between", grp1, "and",
                   grp2, "( n =", length(unique(nodiffTargets$geneId)))
commonGoScatter <- topGO_scatterPlot(df = commonGo, title = plotTitle)

##################################################################################

aligned_plots <- align_plots(
  plotlist = list(tf1 = tf1GoScatter, tf2 = tf2GoScatter, common = commonGoScatter),
  align = "v")

pdf(file = paste(outPrefix, "scatterplot.pdf", sep = ""),
    width = 15, height = 14, onefile = T, pointsize = 18)
ggdraw(aligned_plots$tf1)
ggdraw(aligned_plots$tf2)
ggdraw(aligned_plots$common)
dev.off()


# wd <- (min(max(nchar(as.character(goData$Term))), 80) * 30) * 1.5
wd <- 3500

ht <- max(nrow(tf1Go) * 80, 1500)
png(filename = paste(outPrefixTf1, "_scatterplot.png", sep = ""),
    width = wd, height = ht, res = max(min(wd, ht) / 12, 200))
ggdraw(aligned_plots$tf1)
dev.off()


ht <- max(nrow(tf2Go) * 80, 1500)
png(filename = paste(outPrefixTf2, "_scatterplot.png", sep = ""),
    width = wd, height = ht, res = max(min(wd, ht) / 12, 200))
ggdraw(aligned_plots$tf2)
dev.off()

ht <- max(nrow(commonGo) * 80, 1500)
png(filename = paste(outPrefixCommon, "_scatterplot.png", sep = ""),
    width = wd, height = ht, res = max(min(wd, ht) / 12, 200))
ggdraw(aligned_plots$common)
dev.off()






