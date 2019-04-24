library(chipmine)
library(here)
library(BSgenome.Afumigatus.AspGD.Af293)


rm(list = ls())

outDir <- here::here("analysis", "ChIPseq_analysis", "motif_analysis", "combined_results", "motif_scanning")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

outPrefix <- paste(outDir, "/", "fimo_annotated", sep = "")

##################################################################################

file_fimo <- here::here("analysis", "ChIPseq_analysis", "motif_analysis", "combined_results", "motif_scanning",
                        "fimo_motif_scan_500", "fimo.tsv")

file_peakInfo <- here::here("analysis", "ChIPseq_analysis", "peak_targets", "peakwise_diffbind_groups.tab")

file_motifSource <- here::here("analysis", "ChIPseq_analysis", "motif_analysis", "combined_results",
                               "motif_scanning", "motif_source.txt")

##################################################################################

motifGroups <- suppressMessages(readr::read_tsv(file = file_motifSource, comment = "#"))

fimo <- suppressMessages(readr::read_tsv(file = file_fimo, comment = "#")) %>% 
  dplyr::rename(peakId = sequence_name,
                pValue = `p-value`,
                qValue = `q-value`)

peakInfo <- suppressMessages(readr::read_tsv(file = file_peakInfo)) %>% 
  dplyr::select(peakId, peakPval, pvalFiltered, sampleId, categoryDiffbind,
                summitRegion, geneId, peakPosition) %>% 
  tidyr::separate(col = summitRegion, into = c("chr", "summitSeqStart", "summitSeqEnd"),
                  sep = "(\\:|-)", convert = TRUE)

fimoAnn <- fimo %>% 
  dplyr::left_join(y = motifGroups, by = c("motif_id" = "motifId")) %>% 
  dplyr::left_join(y = peakInfo, by = c("peakId")) %>% 
  dplyr::mutate(motifStart = summitSeqStart + start,
                motifEnd = summitSeqStart + stop) %>% 
  dplyr::arrange(chr, summitSeqStart, pValue) %>% 
  dplyr::select(chr, motifStart, motifEnd, motif_id, score, strand, motif_alt_id, motifGroup,
                pValue, qValue, matched_sequence, peakId, peakPval, categoryDiffbind, geneId,
                start, stop) %>% 
  dplyr::distinct()


readr::write_tsv(x = fimoAnn, path = paste(outPrefix, ".tab", sep = ""))

fimoBed <- dplyr::select(fimoAnn, chr, motifStart, motifEnd, motif_id, score, strand) %>% 
  dplyr::distinct()

readr::write_tsv(x = fimoBed, path = paste(outPrefix, ".bed", sep = ""), col_names = FALSE)








