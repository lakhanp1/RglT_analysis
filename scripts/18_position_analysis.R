library(chipmine)
library(ggpubr)

## plot the raw count and frequency of kmer position analysis results

rm(list = ls())

outDir <- here::here("analysis", "ChIPseq_analysis", "motif_analysis", "combined_results")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

analysisName <- "position_analysis_ccg_gcc"
outPrefix <- paste(outDir, "/", analysisName, sep = "")

##################################################################################

kmersToPlot <- c("ccg", "gcc")

motifAnalysis <- c("CREEHA_CONTROL_enriched" = "CREEHA_CONTROL_enriched",
                   "creE_common" = "CONTROL_10MMAA_common",
                   "CREEHA_10MMAA_enriched" = "CREEHA_10MMAA_enriched")

sampleInfo <- purrr::map_dfr(
  .x = motifAnalysis,
  .f = function(x){
    positionAnalysisDir <- dir(path = here::here("analysis", "ChIPseq_analysis", "motif_analysis", x),
                               pattern = "position_analysis_*", full.names = T)
    list(
      motifGroup = x,
      dir = positionAnalysisDir,
      file = paste(positionAnalysisDir, "/position-analysis.tab", sep = "")
    )
  }, .id = "id"
)

##################################################################################

masterData <- NULL

for (sampleIndex in 1:nrow(sampleInfo)) {
  dt <- suppressMessages(readr::read_tsv(file = sampleInfo$file[sampleIndex], comment = ";")) %>% 
    dplyr::rename(seq = `#seq`)
  
  binCols <- grep(pattern = ".5", x = colnames(dt), value = T)
  
  freqData <- data.table::melt(data = as.data.table(dt), measure.vars = binCols,
                               variable.name = "bin", value.name = "count", variable.factor = FALSE) %>% 
    as.data.frame() %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(total = sum(count),
                  meanCount = mean(count),
                  bin = as.numeric(bin)) %>% 
    dplyr::mutate(frequency = count / total,
                  meanFreq = meanCount / total) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(group = sampleInfo$id[sampleIndex])
  
  if(is.null(masterData)){
    masterData <- freqData
  } else{
    masterData <- dplyr::bind_rows(masterData, freqData)
  }
}

masterData$group <- factor(masterData$group, levels = sampleInfo$id)

masterData <- dplyr::filter(masterData, seq %in% kmersToPlot)


freqPt <- ggplot(data = masterData) +
  geom_line(mapping = aes(x = bin, y = frequency, color = id), size = 1.2) +
  geom_line(mapping = aes(x = bin, y = meanFreq, color = id), alpha = 0.5, size = 0.8, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.8, alpha = 0.5) +
  facet_wrap(facets = ~group, nrow = 3, ncol = 1) +
  scale_color_discrete(name = "3mer") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.x = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 16, face = "bold"),
    panel.grid = element_blank()
  )


countPt <- ggplot(data = masterData) +
  geom_line(mapping = aes(x = bin, y = count, color = id), size = 1.2) +
  geom_line(mapping = aes(x = bin, y = meanCount, color = id), size = 0.8, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2, size = 1.2, alpha = 0.8) +
  facet_wrap(facets = ~group, nrow = 3, ncol = 1, scales = "free_y") +
  scale_color_discrete(name = "3mer") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.x = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 16, face = "bold"),
    panel.grid = element_blank()
  )


mergedPt <- ggarrange(
  countPt, freqPt,
  ncol = 2, common.legend = T, legend = "bottom"
)



png(filename = paste(outPrefix, ".lineplot.png", sep = ""), width = 3500, height = 2000, res = 200)
plot(mergedPt)
dev.off()

pdf(file = paste(outPrefix, ".lineplot.pdf", sep = ""), width = 12, height = 8)
plot(mergedPt)
dev.off()










