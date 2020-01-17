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
endocytic_genes <- unique(read_xls(file.path(proj.dir, "LP_17842 G-CUSTOM-214454.xls")))

# this input derives from the getExpDataMatchedStages.R script
# each dataset contains expression data (counts) 
# on matched tumor/normal samples from patients belonging to different disease stage

patients_earlyStage <- readRDS(file.path(proj.dir, "data_matchedStages", "patients_Earlystages.RDS"))
patients_lateStage <- readRDS(file.path(proj.dir, "data_matchedStages", "patients_Latestages.RDS"))
patients_earlyStageIandII <- readRDS(file.path(proj.dir, "data_matchedStages", "patients_StageIandII.RDS"))
patients_earlyStageIIIandIV <- readRDS(file.path(proj.dir, "data_matchedStages", "patients_StagesIIIandIV.RDS"))

# this input derives from the getExpData.R script
# each dataset contains expression data (counts) 
# on all tumor/normal samples from patients belonging to different disease stage

# patients_earlyStage <- readRDS(file.path(proj.dir, "data", "patients_Earlystages.RDS"))
# patients_lateStage <- readRDS(file.path(proj.dir, "data", "patients_Latestages.RDS"))
# patients_earlyStageIandII <- readRDS(file.path(proj.dir, "data", "patients_StageIandII.RDS"))
# patients_earlyStageIIIandIV <- readRDS(file.path(proj.dir, "data", "patients_StagesIIIandIV.RDS"))

Pattern <- grep("patients", names(.GlobalEnv), value = TRUE)
Patients <- do.call("list", mget(Pattern))


flog.debug("Define dataset for the vizualization")

# p = patients_earlyStage
p = patients_lateStage


flog.debug("Prepare heatmap input for endocytic genes")

HeatmapData <- dplyr::mutate_at(p, vars(-c(Symbol)), funs(log2(.+1)))
Av <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, mean)
StDev <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, sd)
HeatmapData_input <- mutate_at(HeatmapData, vars(-c(Symbol)), funs(((.)-Av)/StDev))
HeatmapData_input <- filter(HeatmapData_input, Symbol %in% endocytic_genes[[1]])

Gene_Symbols <- HeatmapData_input$Symbol
HeatmapData_input2 <- as.data.frame(t(HeatmapData_input[,-1]))
colnames(HeatmapData_input2) <- Gene_Symbols

HeatmapData_input2 <- rownames_to_column(HeatmapData_input2, var = "barcodes")

HeatmapData_input2$status <- with(HeatmapData_input2,
                                   ifelse(test = grepl("-01[A|B]-", HeatmapData_input2$barcodes),
                                          yes = "Tumor", no = "Normal"))

HeatmapData_input2 <- HeatmapData_input2[, c(ncol(HeatmapData_input2), 
                                             1:(ncol(HeatmapData_input2)-1))]
HeatmapData_matrix <- as.matrix(HeatmapData_input2[, c(3:ncol(HeatmapData_input2))])

sample.colors = c("gold", "black")
names(sample.colors) = c("Normal", "Tumor")
status_info = data.frame(status = HeatmapData_input2$status)

fontsize = 0.6


flog.debug("Draw heatmap")
  
Heatmap(t(HeatmapData_matrix),
        cluster_columns = TRUE,
        column_names_side = "top",
        column_dend_side = "top",
        column_names_gp = gpar(cex = fontsize, fontfamily = "Arial"),
        row_names_side = "left",
        #row_dend_side = "left",
        show_row_names = FALSE,
        row_names_gp = gpar(cex = fontsize),
        row_dend_width = unit(3, "cm"),
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D",
        name = "z-score",
        # column_title = "COAD and READ datasets - Early Stages - Tumor vs. Normal",
        column_title = "COAD and READ datasets - Late Stages - Tumor vs. Normal",
        column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
        row_title_rot = 0,
        col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
        row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
        top_annotation = HeatmapAnnotation(
          status = status_info$status,
          col = list(status = c("Tumor" = "gold", "Normal" = "black")),
          show_legend = TRUE,
          show_annotation_name = FALSE))


flog.debug("Prepare heatmap input for ESCRT-I genes")

GOI <- c("insert pool of genes")

HeatmapData <- dplyr::mutate_at(p, vars(-c(Symbol)), funs(log2(.+1)))
Av <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, mean)
StDev <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, sd)
HeatmapData_input <- mutate_at(HeatmapData, vars(-c(Symbol)), funs(((.)-Av)/StDev))
HeatmapData_input <- filter(HeatmapData_input, Symbol %in% GOI)

Gene_Symbols <- HeatmapData_input$Symbol
HeatmapData_input2 <- as.data.frame(t(HeatmapData_input[,-1]))
colnames(HeatmapData_input2) <- Gene_Symbols

HeatmapData_input2 <- rownames_to_column(HeatmapData_input2, var = "barcodes")

HeatmapData_input2$status <- with(HeatmapData_input2,
                                  ifelse(test = grepl("-01[A|B]-", HeatmapData_input2$barcodes),
                                         yes = "Tumor", no = "Normal"))

HeatmapData_input2 <- HeatmapData_input2[, c(ncol(HeatmapData_input2), 
                                             1:(ncol(HeatmapData_input2)-1))]
HeatmapData_matrix <- as.matrix(HeatmapData_input2[, c(3:ncol(HeatmapData_input2))])

sample.colors = c("gold", "black")
names(sample.colors) = c("Normal", "Tumor")
status_info = data.frame(status = HeatmapData_input2$status)

fontsize = 0.6


flog.debug("Draw heatmap")

Heatmap(t(HeatmapData_matrix),
        cluster_columns = TRUE,
        column_names_side = "top",
        column_dend_side = "top",
        column_names_gp = gpar(cex = fontsize, fontfamily = "Arial"),
        row_names_side = "left",
        #row_dend_side = "left",
        show_row_names = TRUE,
        row_names_gp = gpar(cex = fontsize),
        row_dend_width = unit(3, "cm"),
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D",
        name = "z-score",
        # column_title = "COAD and READ datasets - Early Stages - Tumor vs. Normal",
        column_title = "COAD and READ datasets - Late Stages - Tumor vs. Normal",
        column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
        row_title_rot = 0,
        col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
        row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
        top_annotation = HeatmapAnnotation(
          status = status_info$status,
          col = list(status = c("Tumor" = "gold", "Normal" = "black")),
          show_legend = TRUE,
          show_annotation_name = FALSE))
