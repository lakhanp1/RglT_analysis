library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(here)
library(org.AFumigatus.Af293.eg.db)

##
## This script plots the 
## 1) plot1: log2(fold_change) heatmap of specifc genes of interest
## 2) plot2: rlog transformed normalized gene counts heatmap 
## 3) plot1 + plot2 + annotations
##


rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/RNAseq_scripts/DESeq2_functions.R")

analysisName <- "multi_deg"
outDir <- here::here("analysis", "RNAseq_data", "combined_RNAseq")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

outPrefix <- paste(outDir, analysisName, sep = "/")

file_sampleInfo <- here::here("analysis", "RNAseq_data", "sampleInfo.txt")
file_geneInfo <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_version_s03-m05-r09_geneInfo.tab"

degResults <-  c("CEA17_AA_vs_CEA17_C", "5A9_AA_vs_5A9_C",
                 "5A9_C_vs_CEA17_C", "5A9_AA_vs_CEA17_AA")

samples <- c()
plotTitle <- "all DEG comparison"

orgDb <- org.AFumigatus.Af293.eg.db

FDR_cut <- 0.05
lfc_cut <- 0.585
up_cut <- lfc_cut
down_cut <- lfc_cut * -1

####################################################################

diffFiles <- purrr::map_dfr(
  .x = degResults,
  .f = function(x){
    list(
      diffPair = x,
      file_diff = here::here("analysis", "RNAseq_data", x, paste(x, ".DESeq2.tab", sep = "")),
      file_rld = here::here("analysis", "RNAseq_data", x, paste(x, ".rlogCounts.tab", sep = ""))
    )
  })

lfcCol <- "log2FoldChange"
# geneSym <- readr::read_tsv(file = file_geneInfo) %>%
#   distinct(geneId, .keep_all = T)

## or use org.db
geneSym <- AnnotationDbi::select(x = orgDb,
                                 keys = keys(orgDb, keytype = "GID"),
                                 columns = c("DESCRIPTION"),
                                 keytype = "GID") %>% 
  dplyr::rename(geneId = GID)

####################################################################
## import data

## function to extract the log2FoldChange, padj and diff coulumns for each DEG result file
get_foldchange <- function(degFile, name){
  
  degs <- fread(file = degFile, sep = "\t", header = T, stringsAsFactors = F)
  
  newColName <- structure(c(lfcCol, "padj"),
                          names = paste(c("log2.", "padj." ), name, sep = ""))
  
  
  df <- degs %>%
    dplyr::mutate(!! lfcCol := if_else(condition = padj < FDR_cut, true = !! as.name(lfcCol), false = 0)) %>% 
    tidyr::replace_na(purrr::set_names(list(0), nm = c(lfcCol))) %>% 
    dplyr::select(geneId, !! lfcCol, padj) %>%
    dplyr::rename(!!!newColName )
  
  return(df)
}


i <- 1

for(i in 1:nrow(diffFiles)){
  dt <- get_foldchange(degFile = diffFiles$file_diff[i], name = diffFiles$diffPair[i])
  geneSym <- dplyr::left_join(geneSym, dt, by = c("geneId" = "geneId"))
}


# dplyr::filter_at(.tbl = geneSym, .vars = vars(starts_with("padj.")), .vars_predicate = all_vars(is.na(.)))
allData <- dplyr::filter_at(.tbl = geneSym,
                            .vars = vars(starts_with("padj.")),
                            .vars_predicate = any_vars(. < FDR_cut))

####################################################################
rownameCol <- "geneId"
showRowNames <- FALSE
rowNameFontSize <- 14
colNameFontSize <- 14

## fold change heatmap
## fold change heatmap
foldChangeDf <- dplyr::select(allData, !! rownameCol, starts_with("log2")) %>%
  tibble::column_to_rownames(var = rownameCol)

foldChangeMat <- data.matrix(foldChangeDf)

colnames(foldChangeMat) <- gsub(pattern = "log2\\.(.*)", replacement = "\\1", colnames(foldChangeMat), perl = TRUE)


fcHeatmap <- Heatmap(matrix = foldChangeMat,
                     col = colorRamp2(breaks = -3:3, colors = RColorBrewer::brewer.pal(n = 7, name = "PuOr")),
                     cluster_rows = TRUE,
                     clustering_distance_rows = "euclidean",
                     cluster_columns = FALSE,
                     show_row_names = FALSE,
                     row_names_gp = gpar(fontsize = rowNameFontSize),
                     column_names_gp = gpar(fontsize = colNameFontSize), 
                     width = unit(10, "cm"),
                     heatmap_legend_param = list(title = "\nlog2(fold_change)")
)




####################################################################

## log2(fold_change) heatmap with annotation
htList1 <- fcHeatmap

# png(filename = paste(outPrefix, "_fc_heatmap.png", sep = ""), width=4000, height=6000, res = 550)

pdf(file = paste(outPrefix, "_fc_heatmap.pdf", sep = ""), width = 10, height = 10, onefile = TRUE)

draw(object = htList1,
     column_title = plotTitle,
     row_title = "Genes",
     column_title_gp = gpar(fontsize = 14)
)

dev.off()








