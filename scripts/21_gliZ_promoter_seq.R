library(regioneR)
library(org.AFumigatus.Af293.eg.db)
library(TxDb.Afumigatus.Af293.AspGD.GFF)
library(BSgenome.Afumigatus.AspGD.Af293)
library(Biostrings)



rm(list = ls())

analysisName <- "04_gliZ_cluster_promoters"

outDir <- here::here("analysis", "02_ChIPseq_analysis", "05_motif_analysis", analysisName)

outPrefix <- paste(outDir, "/", "gliZ_cluster_promoters", sep = "")


##################################################################################

txdb <- TxDb.Afumigatus.Af293.AspGD.GFF
genome <- BSgenome.Afumigatus.AspGD.Af293
orgDb <- org.AFumigatus.Af293.eg.db
##################################################################################

cdsGrl <- cdsBy(x = txdb, by = "gene")

cdsGr <- unlist(range(cdsGrl))

cdsGr$geneId <- names(cdsGr)

upstreamRegion <- GenomicRanges::promoters(
  x = cdsGr, upstream = 1000, downstream = 0, use.names = TRUE
)

start(upstreamRegion)[which(start(upstreamRegion) < 0)] <- 1
start(upstreamRegion)[which(start(upstreamRegion) == 1)]

upstreamSeq <- BSgenome::getSeq(x = genome, names = upstreamRegion)

upstreamSeq <- upstreamSeq[which(width(upstreamSeq) == 1000)]

gliGenes <- c("Afu6g09630")

gligRegion <- regioneR::toGRanges(A = "Chr6_A_fumigatus_Af293:2,345,058-2,375,633")

glizOvlp <- findOverlaps(query = upstreamRegion, subject = gligRegion)

glizGenes <- upstreamRegion$geneId[glizOvlp@from]

glizGeneProSeq <- upstreamSeq[which(names(upstreamSeq) %in% glizGenes)]

writeXStringSet(x = glizGeneProSeq, filepath = paste(outPrefix, ".fasta", sep = ""))

## random sequences
randomSeq <- upstreamSeq[sample(x = which(!names(upstreamSeq) %in% glizGenes), size = length(glizGenes))]

writeXStringSet(x = randomSeq, filepath = paste(outDir, "/random_promoters.fasta", sep = ""))


AnnotationDbi::select(x = orgDb, keys = glizGenes, columns = "GENE_NAME", keytype = "GID")






