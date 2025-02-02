library(chipmine)
library(here)
library(BSgenome.Afumigatus.Af293.AspGD)
library(data.table)
library(org.AFumigatus.Af293.eg.db)
library(TxDb.Afumigatus.Af293.AspGD.GFF)

## this script
## generates peak target matrix from diffbind annotation results
rm(list = ls())

outDir <- here::here("analysis", "02_ChIPseq_analysis", "01_peak_targets")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

outPrefix <- paste(outDir, "/", sep = "")

compare <- c("CREEHA_CONTROL", "CREEHA_10MMAA")

##################################################################################
summitSeqLen <- 500

file_diffbindInfo <- here::here("analysis", "02_ChIPseq_analysis", "03_diffBind", "sampleInfo.txt")
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_diffbindRes <- here::here("analysis", "02_ChIPseq_analysis", "03_diffBind", "creE_diffbind.all.annotation.tab")
file_diffbindAnn <- here::here("analysis", "02_ChIPseq_analysis", "03_diffBind", "creE_diffbind.annotation.filtered.tab")
TF_dataPath <- here::here("data", "TF_data")

sampleInfo <- suppressMessages(readr::read_tsv(file = file_diffbindInfo))

orgDb <- org.AFumigatus.Af293.eg.db
txDb <- TxDb.Afumigatus.Af293.AspGD.GFF
genomeDb <- BSgenome.Afumigatus.Af293.AspGD

## get the sample details
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   dataPath = TF_dataPath,
                                   profileMatrixSuffix = "normalizedmatrix")

exptDataList <- purrr::transpose(exptData) %>%
  purrr::set_names(nm = purrr::map(., "sampleId"))

grp1 <- compare[1]
grp2 <- compare[2]

grp1Enrich = paste(grp1, ":enriched", sep = "")
grp2Enrich = paste(grp2, ":enriched", sep = "")
grp1Specific = paste(grp1, ":specific", sep = "")
grp2Specific = paste(grp2, ":specific", sep = "")

## prepare ChIPseq target data
bestGrp1Id <- exptData$sampleId[exptData$groupId == grp1 & exptData$bestRep == 1]
bestGrp2Id <- exptData$sampleId[exptData$groupId == grp2 & exptData$bestRep == 1]


tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered", "summitSeq", "summitRegion"),
  FUN = function(x){ structure(paste(x, ".", exptData$sampleId, sep = ""), names = exptData$sampleId) },
  simplify = F, USE.NAMES = T
)

geneStrands <- as.data.frame(GenomicFeatures::genes(x = txDb), row.names = NULL) %>% 
  dplyr::select(geneId = gene_id, strand)

geneInfo <- AnnotationDbi::select(x = orgDb, keys = keys(orgDb), columns = "DESCRIPTION") %>% 
  dplyr::rename(geneId = GID) %>% 
  dplyr::left_join(geneStrands, by = "geneId") %>% 
  dplyr::select(geneId, strand, DESCRIPTION)



##################################################################################

diffbindAnn <- suppressMessages(readr::read_tsv(file = file_diffbindAnn, col_names = T)) %>% 
  dplyr::mutate(DiffBind_region = paste(seqnames, ":", start, "-", end, sep = "")) %>% 
  dplyr::mutate(
    peakAnnotation = if_else(
      condition = is.na(!!sym(tfCols$peakType[bestGrp1Id])),
      true = !!sym(tfCols$peakType[bestGrp2Id]),
      false = !!sym(tfCols$peakType[bestGrp1Id])
    )
  ) %>% 
  dplyr::rename(DiffBind.foldChange = Fold,
                DiffBind.pvalue = p.value,
                DiffBind.FDR = FDR)

peakStats <- dplyr::filter(diffbindAnn, pvalGood.all > 0) %>% 
  dplyr::mutate(peakGene = paste(name, geneId, sep = "::")) %>% 
  dplyr::group_by(categoryDiffbind) %>% 
  dplyr::summarise(count = n_distinct(name))


peakStats$categoryDiffbind <- factor(
  x = peakStats$categoryDiffbind,
  levels =  rev(c(grp1Specific, grp1Enrich, "common", grp2Enrich, grp2Specific))
)

statsPt <- ggplot(data = peakStats, mapping = aes(x = categoryDiffbind, y = count)) +
  geom_bar(stat = "identity", fill = "black", width = 0.7) +
  geom_text(mapping = aes(label = count),
            size = 16, hjust = -0.2, fontface = "bold") +
  scale_y_continuous(expand = expand_scale(add = c(0, 50))) +
  coord_flip() +
  labs(title = "Peak statists") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(size = 2),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_blank(),
    # axis.ticks = element_blank(),
    axis.ticks = element_line(size = 2, linetype = 'dashed'),
    axis.ticks.length = unit(0.3, "cm"),
    axis.title = element_blank(),
    plot.title = element_text(size = 24, face = "bold")
  )

png(filename = paste(outPrefix, "peak_stats.png", sep = ""), width = 1000, height = 1000, res = 80)
statsPt
dev.off()

pdf(file = paste(outPrefix, "peak_stats.pdf", sep = ""), width = 10, height = 10)
statsPt
dev.off()


dplyr::filter(diffbindAnn, pvalGood.all > 0) %>% 
  dplyr::select(geneId, name, peakAnnotation, categoryDiffbind) %>% 
  dplyr::group_by(peakAnnotation) %>% 
  dplyr::summarise(nGenes = n_distinct(geneId),
                   nSites = n_distinct(name))

##################################################################################

tf1SummitSeq <- get_peak_summit_seq(file = exptDataList[[bestGrp1Id]]$peakFile,
                                    peakFormat = "narrowPeak",
                                    sampleId = exptDataList[[bestGrp1Id]]$sampleId,
                                    genome = genomeDb,
                                    length = summitSeqLen)

tf2SummitSeq <- get_peak_summit_seq(file = exptDataList[[bestGrp2Id]]$peakFile,
                                    peakFormat = "narrowPeak",
                                    sampleId = exptDataList[[bestGrp2Id]]$sampleId,
                                    genome = genomeDb,
                                    length = summitSeqLen)

targetMat <- diffbindAnn %>% 
  dplyr::left_join(y = tf1SummitSeq, by = unname(tfCols$peakId[bestGrp1Id])) %>%
  dplyr::left_join(y = tf2SummitSeq, by = unname(tfCols$peakId[bestGrp2Id])) %>%
  dplyr::select(geneId, name, DiffBind_region, DiffBind.foldChange, DiffBind.pvalue, DiffBind.FDR,
                peakAnnotation, peakPosition, diffBind, peakOccupancy, categoryDiffbind, 
                starts_with("peakCall."), starts_with("pvalGood."),
                starts_with("hasPeak."), starts_with("peakPval."), starts_with("peakId."),
                starts_with("summitSeq."), starts_with("summitRegion."),
                starts_with("peakType."), starts_with("peakDist.")) %>% 
  dplyr::select(geneId, name, DiffBind_region, DiffBind.foldChange, DiffBind.pvalue, DiffBind.FDR,
                peakAnnotation, peakPosition, diffBind, peakOccupancy, categoryDiffbind, starts_with("peakCall."),
                contains(bestGrp1Id), contains(bestGrp2Id), starts_with("pvalGood."))

targetMat <- dplyr::left_join(x = targetMat, y = geneInfo, by = "geneId")

readr::write_tsv(x = targetMat, path = paste(outPrefix, "diffbind_allPeak_targets.tab", sep = ""))

##################################################################################
## get individual peak target gene
peakAnnotation <- suppressMessages(readr::read_tsv(file = file_diffbindAnn, col_names = T)) %>% 
  dplyr::select(starts_with("peakId."), geneId, peakPosition) %>% 
  dplyr::filter(!is.na(geneId)) %>% 
  as.data.table() %>% 
  data.table::melt.data.table(
    measure.vars = purrr::map(.x = tfCols, .f = function(x) x[c(bestGrp1Id, bestGrp2Id)]) %>% 
      .[c("peakId")],
    variable.name = "sampleId", value.name = "peakId"
  ) %>% 
  tibble::as_tibble() %>% 
  dplyr::filter(!is.na(peakId)) %>% 
  dplyr::select(-sampleId)

## get diffbind results
diffbindRes <- suppressMessages(readr::read_tsv(file = file_diffbindAnn, col_names = T)) %>% 
  dplyr::select(name, diffBind, peakOccupancy, categoryDiffbind, starts_with("pvalGood."),
                starts_with("peakId."), starts_with("peakPval."), starts_with("pvalFiltered.")) %>% 
  dplyr::left_join(y = tf1SummitSeq, by = structure("name", names = tfCols$peakId[bestGrp1Id])) %>%
  dplyr::left_join(y = tf2SummitSeq, by = structure("name", names = tfCols$peakId[bestGrp2Id])) %>% 
  dplyr::distinct()

diffbindRes <- diffbindRes %>% 
  dplyr::mutate(
    !! sym(tfCols$pvalFiltered[bestGrp1Id]) := if_else(
      condition = !! sym(tfCols$peakPval[bestGrp1Id]) >= exptDataList[[bestGrp1Id]]$pval_cutoff,
      true = TRUE, false = FALSE
    ),
    !! sym(tfCols$pvalFiltered[bestGrp2Id]) := if_else(
      condition = !! sym(tfCols$peakPval[bestGrp2Id]) >= exptDataList[[bestGrp2Id]]$pval_cutoff,
      true = TRUE, false = FALSE
    )
  )

longTargetDf <- data.table::melt.data.table(
  data = as.data.table(diffbindRes),
  measure.vars = purrr::map(.x = tfCols, .f = function(x) x[c(bestGrp1Id, bestGrp2Id)]) %>% 
    .[c("peakId", "peakPval", "pvalFiltered", "summitRegion", "summitSeq")],
  variable.name = "sampleId"
)

levels(longTargetDf$sampleId) <- c(bestGrp1Id, bestGrp2Id)

peakwiseDiff <- as_tibble(longTargetDf) %>% 
  dplyr::filter(!is.na(peakId)) %>% 
  dplyr::left_join(y = peakAnnotation, by = "peakId") %>% 
  dplyr::left_join(y = geneStrands, by = "geneId") %>% 
  dplyr::select(peakId, peakPval, pvalFiltered, sampleId, everything()) %>% 
  dplyr::distinct()

summitSeq <- DNAStringSet(x = peakwiseDiff$summitSeq)
negStrandTargets <- which(peakwiseDiff$strand == "-")
summitSeq[negStrandTargets] <- Biostrings::reverseComplement(summitSeq[negStrandTargets])
peakwiseDiff$summitSeqStranded <- as.character(summitSeq, use.names = FALSE)

readr::write_tsv(x = peakwiseDiff, path = paste(outPrefix, "peakwise_diffbind_groups.tab", sep = ""))


motifData <- dplyr::group_by(peakwiseDiff, name) %>% 
  dplyr::arrange(desc(peakPval)) %>% 
  dplyr::slice(1L) %>% 
  dplyr::ungroup()

readr::write_tsv(x = motifData, path = paste(outPrefix, "summit_seq_data.tab", sep = ""))


