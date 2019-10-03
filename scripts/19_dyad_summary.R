library(chipmine)
library(ggpubr)

## plot the dyad pair preference

rm(list = ls())

outDir <- here::here("analysis", "02_ChIPseq_analysis", "05_motif_analysis", "02_combined_results")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

analysisName <- "dyad_summary"
outPrefix <- paste(outDir, "/", analysisName, sep = "")

##################################################################################

motifAnalysis <- c("CONTROL_enriched" = "02_CREEHA_CONTROL_enriched",
                   "common" = "04_CONTROL_10MMAA_common",
                   # "combined" = "01_creE_combined",
                   "AA_enriched" = "02_CREEHA_10MMAA_enriched")

sampleInfo <- purrr::map_dfr(
  .x = motifAnalysis,
  .f = function(x){
    positionAnalysisDir <- here::here("analysis", "02_ChIPseq_analysis", "05_motif_analysis", "01_meme", x)
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

## add GCC/CCG dyad group
masterData <- dplyr::mutate(
  masterData,
  GC_dyad = dplyr::case_when(
    grepl(pattern = "gcc", x = pair) & grepl(pattern = "ccg", x = pair) ~ "GCC--CCG",
    grepl(pattern = "gcc", x = pair) & !grepl(pattern = "ccg", x = pair) ~ "GCC--NNN",
    grepl(pattern = "ccg", x = pair) & !grepl(pattern = "gcc", x = pair) ~ "CCG--NNN",
    TRUE ~ "NNN--NNN"
  )
)

summary <- dplyr::group_by(masterData, group, GC_dyad) %>% 
  dplyr::summarise(
    count = n(),
    occ_sum = sum(occ)
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    groupDesc = dplyr::case_when(
      group == "AA_enriched" ~ "AA specific\nand\nAA enriched",
      group == "common" ~ "Common",
      group == "CONTROL_enriched" ~ "Control specific\nand\nControl enriched"
    )
  )

summary$group <- factor(summary$group, levels = sampleInfo$id)
summary$GC_dyad <- factor(summary$GC_dyad, levels = c("GCC--NNN", "GCC--CCG", "CCG--NNN", "NNN--NNN"))

pt1 <- ggplot(data = summary) +
  geom_bar(mapping = aes(x = groupDesc, y = occ_sum, fill = GC_dyad),
           stat = "identity", position = position_fill(reverse = TRUE)) +
  labs(y = "dyad pair occurrence in RglT summit sequences") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(
    name = "3mer dyad",
    values = c("GCC--NNN" = "#377eb8", "GCC--CCG" = "#984ea3", "CCG--NNN" = "#4daf4a", "NNN--NNN" = "#999999"),
    breaks = c("GCC--NNN", "GCC--CCG", "CCG--NNN", "NNN--NNN")
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20, hjust = 0.5),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 16, face = "bold"),
    panel.grid = element_blank()
  )

png(filename = paste(outPrefix, ".count_barplot.png", sep = ""), width = 2000, height = 600, res = 150)
plot(pt1)
dev.off()








