library(chipmine)
library(ggpubr)

## plot the dyad pair preference

rm(list = ls())

outDir <- here::here("analysis", "ChIPseq_analysis", "motif_analysis", "combined_results")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

analysisName <- "dyad_summary"
outPrefix <- paste(outDir, "/", analysisName, sep = "")

##################################################################################

motifAnalysis <- c("CREEHA_CONTROL_enriched" = "CREEHA_CONTROL_enriched",
                   "creE_common" = "CONTROL_10MMAA_common",
                   "CREEHA_10MMAA_enriched" = "CREEHA_10MMAA_enriched")

sampleInfo <- purrr::map_dfr(
  .x = motifAnalysis,
  .f = function(x){
    positionAnalysisDir <- here::here("analysis", "ChIPseq_analysis", "motif_analysis", x)
    list(
      motifGroup = x,
      file = paste(positionAnalysisDir, "/dyad_analysis.txt", sep = "")
    )
  }, .id = "id"
)

##################################################################################


masterData <- NULL
sampleIndex <- 1

for (sampleIndex in 1:nrow(sampleInfo)) {
  dt <- suppressMessages(readr::read_tsv(file = sampleInfo$file[sampleIndex], comment = ";")) %>% 
    dplyr::rename(dyad = `#sequence`,
                  id = identifier) %>% 
    dplyr::mutate(pair = gsub(pattern = "n\\{\\d+\\}", replacement = "-", x = id),
                  group = sampleInfo$id[sampleIndex]
    ) 
  
  
  if(is.null(masterData)){
    masterData <- dt
  } else{
    masterData <- dplyr::bind_rows(masterData, dt)
  }
}


summary <- dplyr::group_by(masterData, group) %>% 
  dplyr::summarise(
    gcc = length(which(grepl(pattern = "gcc", x = pair) & !grepl(pattern = "ccg", x = pair))),
    gcc_AND_ccg = length(which(grepl(pattern = "gcc", x = pair) & grepl(pattern = "ccg", x = pair))),
    ccg = length(which(grepl(pattern = "ccg", x = pair) & !grepl(pattern = "gcc", x = pair))),
    others = length(which(!grepl(pattern = "ccg|gcc", x = pair))),
    total = n()
  )


pltData <- data.table::melt(data = as.data.table(summary), measure.vars = c("gcc", "gcc_AND_ccg", "ccg", "others"),
                 variable.name = "kmer", value.name = "dyadCount") %>% 
  as.data.frame()

pltData$group <- factor(pltData$group, levels = sampleInfo$id)

pt <- ggplot(data = pltData) +
  geom_bar(mapping = aes(x = group, y = dyadCount, fill = kmer), stat="identity",
           position = position_fill(reverse = TRUE)) +
  ylab("# dyads") +
  coord_flip() +
  scale_fill_manual(
    name = "3mer dyad",
    values = c("gcc" = "#377eb8", "gcc_AND_ccg" = "#984ea3", "ccg" = "#4daf4a", "others" = "#999999"),
    breaks = c("gcc", "gcc_AND_ccg", "ccg", "others"),
    labels = c("GCC--NNN", "GCC--CCG", "CCG--NNN", "NNN--NNN")
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 16, face = "bold"),
    panel.grid = element_blank()
  )

png(filename = paste(outPrefix, ".count_barplot.png", sep = ""), width = 2000, height = 500, res = 150)
plot(pt)
dev.off()








