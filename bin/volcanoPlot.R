rm(list = ls())

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("TCGAbiolinks")

suppressPackageStartupMessages({
  library(futile.logger)
  library(readxl)
  library(WriteXLS)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
})

flog.threshold(DEBUG)

flog.debug("Set working directory and load required data")
proj.dir = "~/Desktop/HT projects/TCGA_transcriptomics_proj"
spec.dir = "data_matchedStages"
endocytic_genes <- unique(read_xls(file.path(proj.dir, "LP_17842 G-CUSTOM-214454.xls")))

file.list <- list.files(file.path(proj.dir, spec.dir), pattern = "dataDEG_")

for (f in file.list) {
  temp.file.name = str_sub(f, 9, str_length(f)-4)
  assign(paste0("DGE_", temp.file.name),
         readRDS(file.path(proj.dir, spec.dir, f)))
}


Pattern <- grep("DGE_", names(.GlobalEnv), value = TRUE)
list.dataDEGs <- do.call("list", mget(Pattern))

GOI <- c("insert pool of genes")

dataDEGs <- DGE_Earlystages
# dataDEGs <- DGE_Latestages

GOI_filtered <- filter(dataDEGs, dataDEGs$Symbol %in% GOI)
vp <- ggplot(data = dataDEGs,
             aes(x = logFC,
             y = -log10(FDR),
             colour = Threshold)) +
  scale_color_manual(values = c("dodgerblue", "gold", "deeppink2")) +
  geom_point(alpha = 0.4, size = 1.0) + xlim(c(-3.5, 3.5)) + ylim(c(0, 10)) +
  labs(color = "Expression pattern") +
  geom_point(data = GOI_filtered, colour = "black") +
  ggtitle(paste0("COAD and READ datasets - Early stages - Tumor vs. Normal")) +
  theme(plot.title = element_text(face = "bold")) +
  xlab("log2FoldChange") + ylab("-log10(FDR)") +
  theme(plot.title = element_text(hjust = 0.5))

vp + geom_text_repel(data = GOI_filtered,
                     aes(label = GOI_filtered$Symbol), colour = "black", size = 4)