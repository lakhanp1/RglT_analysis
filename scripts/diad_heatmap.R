library(chipmine)
library(org.AFumigatus293.eg.db)
library(seriation)


rm(list = ls())

outDir <- here::here("analysis", "ChIPseq_analysis", "motifAnalysis", "rsat_diad_creE_CONTROL_10MMAA_combined")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

outPrefix <- paste(outDir, "/diad_motif_occurrence.all200", sep = "")
##################################################################################
orgDb <- org.AFumigatus293.eg.db

file_diadMatch <- here::here("analysis", "ChIPseq_analysis", "motifAnalysis", "rsat_diad_creE_CONTROL_10MMAA_combined",
                             "creE_CONTROL_10MMAA_combined.diad_match.tab")

file_peakGroup <- here::here("analysis", "ChIPseq_analysis", "motifAnalysis", "rsat_diad_creE_CONTROL_10MMAA_combined",
                             "peak_group.txt")

file_diad <- here::here("analysis", "ChIPseq_analysis", "motifAnalysis", "rsat_diad_creE_CONTROL_10MMAA_combined",
                        "diad_summary.txt")

##################################################################################

peakGroups <- readr::read_tsv(file = file_peakGroup)
diadSummary <- suppressMessages(readr::read_tsv(file = file_diad)) %>% 
  dplyr::mutate(diad = gsub(pattern = "(n\\{|\\})", replacement = "_", x = sequence))

diadMap <- data.table::fread(file = file_diadMatch, sep = "\t") %>% 
  dplyr::filter(PatID != "START_END") %>% 
  # dplyr::filter(Start >= 50, End <= 150) %>%
  # dplyr::filter(Start < 50 | End > 150) %>%
  dplyr::mutate(
    relativeStart = Start / 200,
    diad = gsub(pattern = "(n\\{|\\})", replacement = "_", x = Pattern),
    condition = gsub(pattern = "_withCtrl_peak_\\d+", replacement = "", x = SeqID, perl = T)
  ) %>% 
  dplyr::mutate(summitDist = 0.5 - relativeStart) %>% 
  dplyr::mutate(absSummitDist = abs(summitDist)) %>% 
  dplyr::mutate(invSummitDist = 0.5 - absSummitDist)


groupStats <- diadMap %>% 
  dplyr::filter(Start >= 50, End <= 150) %>%
  dplyr::group_by(condition, Pattern) %>% 
  dplyr::summarise(n = n()) %>% 
  tidyr::spread(condition, n) %>% 
  dplyr::mutate(
    total_10MMAA4 = sum(CREEHA_10MMAA4),
    total_CONTROL4 = sum(CREEHA_CONTROL4)
  )

diadPosDf <- dplyr::group_by(.data = diadMap, SeqID, diad) %>% 
  dplyr::arrange(absSummitDist, .by_group = TRUE) %>% 
  dplyr::slice(1L) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(SeqID, diad, invSummitDist) %>% 
  tidyr::spread(key = diad, value = invSummitDist, fill = 0) %>% 
  dplyr::left_join(y = peakGroups, by = c("SeqID" = "peakId")) %>% 
  dplyr::filter(!is.na(diffbind))


diadPosMat <- dplyr::select(diadPosDf, -SeqID, -diffbind, -group) %>% 
  as.matrix()

rowAnDf <- dplyr::select(diadPosDf, SeqID, group, diffbind) %>% 
  as.data.frame()

rowAnDf$group <- factor(rowAnDf$group, levels = sort(unique(rowAnDf$group)))
rowAnDf$diffbind <- factor(rowAnDf$diffbind, levels = sort(unique(rowAnDf$diffbind)))

rownames(diadPosMat) <- diadPosDf$SeqID

colAnDf <- data.frame(diad = colnames(diadPosMat), stringsAsFactors = F) %>% 
  dplyr::left_join(y = diadSummary, by = "diad") %>% 
  dplyr::select(diad, pair, rank, gap, up, down)


seriateOrd = seriate(max(diadPosMat) - diadPosMat, method = "BEA_TSP")

o1 = seriate(dist(diadPosMat), method = "TSP")
o2 = seriate(dist(t(diadPosMat)), method = "TSP")

mainHt <- Heatmap(matrix = diadPosMat,
                  name = "main",
                  col = colorRamp2(breaks = c(0, 0.5), colors = c("black", "green")),
                  show_row_names = F,
                  row_order = get_order(o1),
                  column_order = get_order(o2),
                  row_split = rowAnDf$diffbind, cluster_row_slices = FALSE,
                  column_split = colAnDf$gap, cluster_column_slices = FALSE,
                  row_title = NULL,
                  column_title = NULL)

groupAn <- rowAnnotation(
  df = rowAnDf[, "diffbind", drop = F],
  name = "diffbind",
  col = list(diffbind = structure(RColorBrewer::brewer.pal(5, "Dark2"), names = levels(rowAnDf$diffbind)))
)

htList <- groupAn + mainHt

pdf(file = paste(outPrefix, ".gap.heatmap.pdf", sep = ""), width = 16, height = 10)

draw(htList,
     main_heatmap = "main",
     column_title = "diad pairs in motif sequences grouped by diffbind: split columns by gap"
)
dev.off()





