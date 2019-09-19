library(chipmine)
library(org.AFumigatus.Af293.eg.db)
library(TxDb.Afumigatus.Af293.AspGD.GFF)
library(here)
library(ggbeeswarm)
library(ggpubr)
library(ggrepel)


rm(list = ls())

##################################################################################
analysisName <- "TF_replicate"
outDir <- here::here("analysis", "02_ChIPseq_analysis", "01_correlation_analysis")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_replicates <- here::here("analysis", "02_ChIPseq_analysis", "01_correlation_analysis", "tf_replicates.txt")

file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")

orgDb <- org.AFumigatus.Af293.eg.db
txDb <- TxDb.Afumigatus.Af293.AspGD.GFF

TF_dataPath <- here::here("data", "TF_data")

file_tf_macs2 <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")

matrixType <- "2kb_summit"
up <- 2000
down <- 2000
body <- 0
bin <- 10
matrixDim = c(c(up, body, down)/bin, bin)

##################################################################################

tfSampleList <- readr::read_tsv(file = file_tf_macs2, col_names = c("id"),  comment = "#")

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSampleList$id,
  dataPath = TF_dataPath,
  profileMatrixSuffix = matrixType)

# tfInfoList <- purrr::transpose(tfInfo)  %>% 
#   purrr::set_names(nm = purrr::map(., "sampleId"))


replicatePairs <- suppressMessages(readr::read_tsv(file = file_replicates))

plotListAll <- list()
plotList_pval_distibution <- list()
i <- 1

for (i in 1:nrow(replicatePairs)) {
  print(replicatePairs$rep1[i])
  
  repInfo <- dplyr::filter(tfInfo, sampleId %in% c(replicatePairs$rep1[i], replicatePairs$rep2[i]))
  
  plots_pval <- tf_replicate_plots(sampleInfo = repInfo, compare = "pvalue",
                                   title = "set1", yintercept = 25)
  
  # plots_enrichment <- tf_replicate_plots(sampleInfo = repInfo, compare = "enrichment",
  #                                        title = "set1", yintercept = 3)
  
  ## all plots combined in a row for a replicate
  repPlot <- ggarrange(
    plots_pval$table, plots_pval$distribution, plots_pval$valueScatter, plots_pval$rankScatter,
    # plots_enrichment$distribution, plots_enrichment$valueScatter, plots_enrichment$rankScatter,
    nrow = 1) +
    theme(plot.background = element_rect(color = "black"),
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))

  plotListAll[[i]] <- repPlot
  
  ## summary table and distribution plot combined
  pvalDistrPlot <- ggarrange(
    plots_pval$table, plots_pval$distribution,
    ncol = 1, heights = c(1, 8)
  ) +
    theme(plot.background = element_rect(color = "black"),
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))
  
  plotList_pval_distibution[[i]] <- pvalDistrPlot
}


figureAll <- ggarrange(plotlist = plotListAll,
                       ncol = 1)


png(filename = paste(outPrefix, ".correlation_pval_enrichment.png", sep = ""),
    width = 7000, height = 3000, res = 280)
figureAll
dev.off()


fig_pval_distribution <- ggarrange(plotlist = plotList_pval_distibution,
                                   nrow = 1, ncol = 2, hjust = 0.5)

png(filename = paste(outPrefix, ".pval_distribution.png", sep = ""),
    width = 4000, height = 2500, res = 350)
fig_pval_distribution
dev.off()






