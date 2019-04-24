library(chipmine)
library(org.AFumigatus293.eg.db)
library(cowplot)

rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/GO_enrichment/topGO_functions.R")

outDir <- here::here("analysis", "ChIPseq_analysis")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

##################################################################################

file_goMap <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_orgDb/geneid2go.AFumigatus_Af293.topGO.map"
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_version_s03-m05-r06_CDS_ORFs.bed"

TF_dataPath <- here::here("data", "TF_data")

sampleList <- c("CREEHA_CONTROL4", "CREEHA_CONTROL5", "CREEHA_10MMAA4", "CREEHA_10MMAA5")
file_targets <- here::here("analysis", "ChIPseq_analysis", "peak_targets", "peak_targets.curated.filtered.tab")
orgDb <- org.AFumigatus293.eg.db
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
        "peakCoverage", "pvalFiltered"),
  FUN = function(x){ structure(paste(x, ".", sampleList, sep = ""), names = sampleList) },
  simplify = F, USE.NAMES = T
)

targets <- suppressMessages(readr::read_tsv(file = file_targets, col_names = T))


##################################################################################


s1 <- "CREEHA_CONTROL4"
s2 <- "CREEHA_10MMAA4"

## TF1 target list
s1Targets <- dplyr::filter_at(.tbl = targets, .vars = tfCols$hasPeak[s1], .vars_predicate = all_vars(. == TRUE))

nrow(s1Targets)
length(unique(s1Targets$geneId))
dplyr::group_by_at(.tbl = s1Targets, .vars = vars(starts_with("hasPeak."), starts_with("pvalFiltered"))) %>% 
  dplyr::summarise(n = length(unique(geneId)))


## TF2 target list
s2Targets <- dplyr::filter_at(.tbl = targets, .vars = tfCols$hasPeak[s2], .vars_predicate = all_vars(. == TRUE))

nrow(s2Targets)
length(unique(s2Targets$geneId))


##################################################################################

## genes specific to CREEHA_CONTROL sample
nodeSize <- 5
plotTitle <- paste(s1, "specific targets")
outPrefix <- paste(outDir, "/topGO/topGO_", s1, "_specific", sep = "")

tf1SpecificTargets <- setdiff(s1Targets$geneId, s2Targets$geneId)

tf1Go <- topGO_enrichment(genes = tf1SpecificTargets, goMapFile = file_goMap, goNodeSize = nodeSize)


tf1Go$condition <- paste(s1, "specific")

readr::write_tsv(x = tf1Go, path = paste(outPrefix, ".tab", sep = ""))

tf1GoScatter <- topGO_scatterPlot(df = tf1Go, title = plotTitle)

plotOut <- paste(outPrefix, ".pdf", sep = "")



##################################################################################
## genes specific to CREEHA_10MMAA3 sample

plotTitle <- paste(s2, "specific targets")
outPrefix <- paste(outDir, "/topGO/topGO_", s2, "_specific", sep = "")

tf2SpecificTargets <- setdiff(s2Targets$geneId, s1Targets$geneId)

tf2Go <- topGO_enrichment(genes = tf2SpecificTargets, goMapFile = file_goMap, goNodeSize = nodeSize)

tf2Go$condition <- paste(s2, "specific")

readr::write_tsv(x = tf2Go, path = paste(outPrefix, ".tab", sep = ""))

tf2GoScatter <- topGO_scatterPlot(df = tf2Go, title = plotTitle)

plotOut <- paste(outPrefix, ".pdf", sep = "")



##################################################################################
## common targets for CREEHA_CONTROL2 and CREEHA_10MMAA3 samples

plotTitle <- paste("common targets between", s1, "and", s2, sep = " ")
outPrefix <- paste(outDir, "/topGO/topGO_", s1, ".", s2, "_common", sep = "")

commonTargets <- intersect(s1Targets$geneId, s2Targets$geneId)

commonGo <- topGO_enrichment(genes = commonTargets, goMapFile = file_goMap, goNodeSize = nodeSize)

commonGo$condition <- "common targets"

readr::write_tsv(x = commonGo, path = paste(outPrefix, ".tab", sep = ""))

commonGoScatter <- topGO_scatterPlot(df = commonGo, title = plotTitle)

plotOut <- paste(outPrefix, ".pdf", sep = "")


##################################################################################

aligned_plots <- align_plots(
  plotlist = list(tf1 = tf1GoScatter, tf2 = tf2GoScatter, common = commonGoScatter),
  align = "v")

pdf(file = paste(outDir, "/topGO/topGO_enrichment.pdf", sep = ""), width = 12, height = 12, onefile = T, pointsize = 18)
ggdraw(aligned_plots$tf1)
ggdraw(aligned_plots$tf2)
ggdraw(aligned_plots$common)
dev.off()



