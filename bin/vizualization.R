################################
# export filtered data frames

# dfFILT.file.cases.all

write(dfFILT.file.cases.all, "dfFILT.file.cases.all.txt", sep = ",")

lst_dfFILT <- mget(dfFILT.file.cases.all)

sapply(1:length(lst_dfFILT),
       function(i) write.table(lst_dfFILT[[i]], file = paste0(dfFILT.file.cases.all[i], ".txt")))

# dfFILT.rep_patient.file.cases.all
lst_dfFILT_rep_patient <- mget(dfFILT.rep_patient.file.cases.all)

heatmap_rep_patients <- reduce(lst_dfFILT_rep_patient, full_join, by = 'Symbol')

# get the file name numbers
prefix <- c("NT", "TP")
suffix <- gsub("dataFilt_rep_patient_TCGA-", "", dfFILT.rep_patient.file.cases.all)
rep_patient_col_names <- paste(prefix, rep(suffix, each = 2), sep = "_")

# substitute colnames
df_temp <- heatmap_rep_patients[, c(2:41)]
colnames(df_temp) <- rep_patient_col_names

# extract order
df_order <- heatmap_rep_patients[, 1]

# combine dataframes
heatmap_rep_patients <- cbind(df_order, df_temp)
heatmap_rep_patients <- dplyr::rename(heatmap_rep_patients, "Symbol" = "df_order")
is.na(heatmap_rep_patients)

write.table(heatmap_rep_patients, "TCGA_master_heatmap_rep_patients.txt", sep = "\t", quote = FALSE)


######################################################################
# read data frame with reads in genes for different cancer types

heatmap_rep_patients <- read.table("TCGA_master_heatmap_rep_patients.txt", header = TRUE)
heatmap_rep_patients[is.na(heatmap_rep_patients)] <- 0

# heatmap of endocytic modulators across human cancer types

endocytic_genes <- read_xls("LP_17842 G-CUSTOM-214454.xls")

HeatmapData <- dplyr::mutate_at(heatmap_rep_patients, vars(-c(Symbol)), funs(log2(.+1)))
Av <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, mean)
StDev <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, sd)
HeatmapData_input <- mutate_at(HeatmapData, vars(-c(Symbol)), funs(((.)-Av)/StDev))
HeatmapData_input <- filter(HeatmapData_input, Symbol %in% endocytic_genes[[5]])

# dim(HeatmapData_input)
# View(HeatmapData_input)

Gene_Symbols <- HeatmapData_input$Symbol
HeatmapData_input2 <- as.data.frame(t(HeatmapData_input[,-1]))
colnames(HeatmapData_input2) <- Gene_Symbols

HeatmapData_input2 <- rownames_to_column(HeatmapData_input2, var = "project")

HeatmapData_input2$status <- with(HeatmapData_input2,
                                  ifelse(test = grepl("TP_", 
                                                      HeatmapData_input2$project),
                                         yes = "Tumor", no = "Normal"))

HeatmapData_input2 <- HeatmapData_input2[, c(ncol(HeatmapData_input2), 
                                             1:(ncol(HeatmapData_input2)-1))]

# remove redundant symbols TP_ and NT_

# temp_input <- str_replace_all(HeatmapData_input2$project, "TP\\_|NT\\_", "")
# HeatmapData_input2$project <- temp_input

# dim(HeatmapData_input2)
# View(HeatmapData_input2)

HeatmapData_matrix <- as.matrix(HeatmapData_input2[, c(3:ncol(HeatmapData_input2))])

sample.colors <- c("gold", "black")
sample.colors 

names(sample.colors) <- c("Normal", "Tumor")
sample.colors

status_info <- data.frame(status = HeatmapData_input2$status)

# View(t(HeatmapData_matrix))


Heatmap(t(HeatmapData_matrix),
        cluster_columns = TRUE,
        column_names_side = "top",
        column_dend_side = "top",
        column_names_gp = gpar(cex = fontsize, fontfamily = "Arial"),
        row_names_side = "left",
        show_row_names = FALSE,
        row_names_gp = gpar(cex = fontsize),
        row_dend_width = unit(3, "cm"),
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D",
        name = "z-score",
        #column_title = "Expression of endocytic genes in across cancer types",
        column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
        row_title_rot = 0,
        col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
        row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
        bottom_annotation = HeatmapAnnotation(status_info,
                                              col = list(status = sample.colors),
                                              show_legend = TRUE))



######################################################################
# read data frame with reads for genes in different cancer types

heatmap_rep_patients <- read.table("TCGA_master_heatmap_rep_patients.txt", header = TRUE)
heatmap_rep_patients[is.na(heatmap_rep_patients)] <- 0

# heatmap of endocytic modulators across human cancer types

endocytic_genes <- read_xls("LP_17842 G-CUSTOM-214454.xls")

HeatmapData <- dplyr::mutate_at(heatmap_rep_patients, vars(-c(Symbol)), funs(log2(.+1)))
Av <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, mean)
StDev <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, sd)
HeatmapData_input <- mutate_at(HeatmapData, vars(-c(Symbol)), funs(((.)-Av)/StDev))
HeatmapData_input <- filter(HeatmapData_input, Symbol %in% endocytic_genes[[5]])

# dim(HeatmapData_input)

Gene_Symbols <- HeatmapData_input$Symbol
HeatmapData_input2 <- as.data.frame(t(HeatmapData_input[,-1]))
colnames(HeatmapData_input2) <- Gene_Symbols
HeatmapData_input2 <- rownames_to_column(HeatmapData_input2, var = "project")

HeatmapData_input2$status <- with(HeatmapData_input2,
                                  ifelse(test = grepl("TP_", 
                                                      HeatmapData_input2$project),
                                         yes = "Tumor", no = "Normal"))

HeatmapData_input2 <- HeatmapData_input2[, c(ncol(HeatmapData_input2), 
                                             1:(ncol(HeatmapData_input2)-1))]

status_info <- data.frame(status = HeatmapData_input2$status)

#temp_input <- gsub("TP_", "", HeatmapData_input2$project)
#temp_input <- gsub("NT_", ".", temp_input)

#HeatmapData_input2$project <- temp_input

HeatmapData_input2 <- column_to_rownames(HeatmapData_input2, var = "project")
HeatmapData_input2 <- HeatmapData_input2[, c(2:ncol(HeatmapData_input2))]


# dim(HeatmapData_input2)

HeatmapData_matrix <- as.matrix(HeatmapData_input2)

sample.colors <- c("gold", "black")
sample.colors 

names(sample.colors) <- c("Normal", "Tumor")
sample.colors

Heatmap(t(HeatmapData_matrix),
        cluster_columns = TRUE,
        column_names_side = "top",
        column_dend_side = "top",
        column_names_gp = gpar(cex = fontsize, fontfamily = "Arial"),
        row_names_side = "left",
        show_row_names = FALSE,
        row_names_gp = gpar(cex = fontsize),
        row_dend_width = unit(3, "cm"),
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D",
        name = "z-score",
        #column_title = "Expression of endocytic genes in across cancer types",
        column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
        row_title_rot = 0,
        col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
        row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
        top_annotation = HeatmapAnnotation(status_info,
                                           col = list(status = sample.colors),
                                           show_legend = TRUE))



#######################################################################################
# read data frame with reads for genes in different cancer types (normal vs. cancer)

heatmap_rep_patients <- read.table("TCGA_master_heatmap_rep_patients.txt", header = TRUE)
heatmap_rep_patients[is.na(heatmap_rep_patients)] <- 0

# heatmap of endocytic modulators across human cancer types

endocytic_genes <- read_xls("LP_17842 G-CUSTOM-214454.xls")

HeatmapData <- dplyr::mutate_at(heatmap_rep_patients, vars(-c(Symbol)), funs(log2(.+1)))
Av <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, mean)
StDev <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, sd)
HeatmapData_input <- mutate_at(HeatmapData, vars(-c(Symbol)), funs(((.)-Av)/StDev))
HeatmapData_input <- filter(HeatmapData_input, Symbol %in% endocytic_genes[[5]])

Gene_Symbols <- HeatmapData_input$Symbol
HeatmapData_input2 <- as.data.frame(t(HeatmapData_input[,-1]))
colnames(HeatmapData_input2) <- Gene_Symbols

HeatmapData_input2 <- rownames_to_column(HeatmapData_input2, var = "project")

HeatmapData_input2$status <- with(HeatmapData_input2,
                                  ifelse(test = grepl("TP_", 
                                                      HeatmapData_input2$project),
                                         yes = "Tumor", no = "Normal"))

HeatmapData_input2 <- HeatmapData_input2[, c(ncol(HeatmapData_input2), 
                                             1:(ncol(HeatmapData_input2)-1))]

#status_info <- data.frame(status = HeatmapData_input2$status)

tumour <- filter(HeatmapData_input2, status == "Tumor")
tumour <- tumour[, -1]
tumour <- column_to_rownames(tumour, var = "project")
tumour <- as.matrix(tumour)
normal <- filter(HeatmapData_input2, status == "Normal")
normal <- normal[, -1]
normal <- column_to_rownames(normal, var = "project")
normal <- as.matrix(normal)

ht1 = Heatmap(t(normal),
              column_title = "Normal",
              cluster_columns = TRUE,
              column_names_side = "top",
              column_dend_side = "top",
              column_names_gp = gpar(cex = fontsize, fontfamily = "Arial"),
              row_names_side = "left",
              show_row_names = FALSE,
              row_dend_side = "left",
              row_names_gp = gpar(cex = fontsize),
              row_dend_width = unit(1.5, "cm"),
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "ward.D",
              name = "z-score",
              column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
              row_title_rot = 0,
              col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")), 
              row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"))

#ht1

ht2 = Heatmap(t(tumour),
              column_title = "Tumor",
              cluster_columns = TRUE,
              column_names_side = "top",
              column_dend_side = "top",
              column_names_gp = gpar(cex = fontsize, fontfamily = "Arial"),
              row_names_side = "left",
              show_row_names = FALSE,
              row_dend_side = "left",
              row_names_gp = gpar(cex = fontsize),
              row_dend_width = unit(1.5, "cm"),
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "ward.D",
              name = "z-score",
              column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
              row_title_rot = 0,
              col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
              row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"))

#ht2

ht1 + ht2


# heatmap of ESCRT genes across human cancer types

endocytic_genes <- c("HGS", "STAM", "STAM2",
                     "TSG101", "VPS28", "VPS37A", "VPS37B", "VPS37C", "VPS37D", 
                     "UBAP1", "MVB12A", "MVB12B",
                     "VPS36", "SNF8", "VPS25",
                     "CHMP6", "CHMP4A", "CHMP4B", "CHMP4C",
                     "CHMP3", "CHMP2A", "CHMP2B", "CHMP5", "PDCD6IP",
                     "IST1", "CHMP1", "VPS4A", "VPS4B")

HeatmapData <- dplyr::mutate_at(heatmap_rep_patients, vars(-c(Symbol)), funs(log2(.+1)))
Av <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, mean)
StDev <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, sd)
HeatmapData_input <- mutate_at(HeatmapData, vars(-c(Symbol)), funs(((.)-Av)/StDev))
HeatmapData_input <- filter(HeatmapData_input, Symbol %in% endocytic_genes)

# dim(HeatmapData_input)

Gene_Symbols <- HeatmapData_input$Symbol
HeatmapData_input2 <- as.data.frame(t(HeatmapData_input[,-1]))
colnames(HeatmapData_input2) <- Gene_Symbols
HeatmapData_input2 <- rownames_to_column(HeatmapData_input2, var = "project")

HeatmapData_input2$status <- with(HeatmapData_input2,
                                  ifelse(test = grepl("TP_", 
                                                      HeatmapData_input2$project),
                                         yes = "Tumor", no = "Normal"))

HeatmapData_input2 <- HeatmapData_input2[, c(ncol(HeatmapData_input2), 
                                             1:(ncol(HeatmapData_input2)-1))]

status_info <- data.frame(status = HeatmapData_input2$status)

HeatmapData_input2 <- column_to_rownames(HeatmapData_input2, var = "project")
HeatmapData_input2 <- HeatmapData_input2[, c(2:ncol(HeatmapData_input2))]

# dim(HeatmapData_input2)

HeatmapData_matrix <- as.matrix(HeatmapData_input2)

sample.colors <- c("gold", "black")
sample.colors 

names(sample.colors) <- c("Normal", "Tumor")
sample.colors

Heatmap(t(HeatmapData_matrix),
        cluster_columns = TRUE,
        column_names_side = "top",
        column_dend_side = "top",
        column_names_gp = gpar(cex = fontsize, fontfamily = "Arial"),
        row_names_side = "left",
        show_row_names = TRUE,
        row_names_gp = gpar(cex = fontsize),
        row_dend_width = unit(3, "cm"),
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D",
        name = "z-score",
        #column_title = "Expression of endocytic genes in across cancer types",
        column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
        row_title_rot = 0,
        col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
        row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
        top_annotation = HeatmapAnnotation(status_info,
                                           col = list(status = sample.colors),
                                           show_legend = TRUE))

# Different repressentation

HeatmapData <- dplyr::mutate_at(heatmap_rep_patients, vars(-c(Symbol)), funs(log2(.+1)))
Av <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, mean)
StDev <- apply(HeatmapData[, c(2:ncol(HeatmapData))], 1, sd)
HeatmapData_input <- mutate_at(HeatmapData, vars(-c(Symbol)), funs(((.)-Av)/StDev))
HeatmapData_input <- filter(HeatmapData_input, Symbol %in% endocytic_genes)

Gene_Symbols <- HeatmapData_input$Symbol
HeatmapData_input2 <- as.data.frame(t(HeatmapData_input[,-1]))
colnames(HeatmapData_input2) <- Gene_Symbols

HeatmapData_input2 <- rownames_to_column(HeatmapData_input2, var = "project")

HeatmapData_input2$status <- with(HeatmapData_input2,
                                  ifelse(test = grepl("TP_", 
                                                      HeatmapData_input2$project),
                                         yes = "Tumor", no = "Normal"))

HeatmapData_input2 <- HeatmapData_input2[, c(ncol(HeatmapData_input2), 
                                             1:(ncol(HeatmapData_input2)-1))]

tumour <- filter(HeatmapData_input2, status == "Tumor")
tumour <- tumour[, -1]
tumour <- column_to_rownames(tumour, var = "project")
tumour <- as.matrix(tumour)
normal <- filter(HeatmapData_input2, status == "Normal")
normal <- normal[, -1]
normal <- column_to_rownames(normal, var = "project")
normal <- as.matrix(normal)

ht1 = Heatmap(t(normal),
              column_title = "Normal",
              cluster_columns = TRUE,
              column_names_side = "top",
              column_dend_side = "top",
              column_names_gp = gpar(cex = fontsize, fontfamily = "Arial"),
              row_names_side = "left",
              show_row_names = TRUE,
              row_dend_side = "left",
              row_names_gp = gpar(cex = fontsize),
              row_dend_width = unit(1.5, "cm"),
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "ward.D",
              name = "z-score",
              column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
              row_title_rot = 0,
              col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")), 
              row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"))

#ht1

ht2 = Heatmap(t(tumour),
              column_title = "Tumor",
              cluster_columns = TRUE,
              column_names_side = "top",
              column_dend_side = "top",
              column_names_gp = gpar(cex = fontsize, fontfamily = "Arial"),
              row_names_side = "left",
              show_row_names = FALSE,
              row_dend_side = "left",
              row_names_gp = gpar(cex = fontsize),
              row_dend_width = unit(1.5, "cm"),
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "ward.D",
              name = "z-score",
              column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"),
              row_title_rot = 0,
              col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
              row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"))

#ht2

ht1 + ht2

################################
# export filtered data frames

# dfFILT.file.cases.all

write(dfDEG.file.cases.all, "dfDEG.file.cases.all.txt", sep = ",")

lst_dfDEG <- mget(dfDEG.file.cases.all)

sapply(1:length(lst_dfDEG),
       function(i) write.table(lst_dfDEG[[i]], file = paste0(dfDEG.file.cases.all[i], ".txt")))


#####################################################
# prepare and export a master data frame with DEGs

#dfDEG.file.cases.all
lst_dfDEG <- mget(dfDEG.file.cases.all)

dfDEG <- reduce(lst_dfDEG, full_join, by = 'Symbol')

# get the file name numbers
prefix <- c("logFC", "logCPM", "LR", "PValue", "FDR")
suffix <- gsub("dfDEG_TCGA-", "", dfDEG.file.cases.all)
DEGs_col_names <- paste(prefix, rep(suffix, each = 5), sep = "_")

# substitute colnames
df_temp <- dfDEG[, c(2:101)]
colnames(df_temp) <- DEGs_col_names

# extract order
df_order <- dfDEG[, 1]

# combine dataframes
dfDEG <- cbind(df_order, df_temp)
dfDEG <- dplyr::rename(dfDEG, "Symbol" = "df_order")
is.na(dfDEG)
dfDEG[is.na(dfDEG)] <- 0
is.na(dfDEG)

write.table(dfDEG, "TCGA_master_dfDEGs.txt", sep = "\t", quote = FALSE)


dfDEG <- read.table("TCGA_master_dfDEGs.txt", sep = "\t", row.names = NULL)

dfDEG_selected <- dplyr::select(dfDEG, "Symbol", starts_with("logFC"), starts_with("FDR"))
View(dfDEG_selected)
dim(dfDEG_selected)

# changes across different cancer types

endocytic_genes <- read_xls("LP_17842 G-CUSTOM-214454.xls")
endocytic_genes <- endocytic_genes$`Gene Symbol` 

dfDEG_selected_endocytic_genes <- filter(dfDEG_selected, Symbol %in% endocytic_genes)

consistently_up <- dfDEG_selected_endocytic_genes[, c(1:20)] %>%
  filter_if(is.numeric, all_vars(. >= 0.323))

alt_up <- dfDEG_selected_endocytic_genes[, c(1:20)] %>%
  filter_if(is.numeric, any_vars(. >= 0.323))

consistently_down <- dfDEG_selected_endocytic_genes[, c(1:20)] %>%
  filter_if(is.numeric, all_vars(. <= -0.323))

alt_down <- dfDEG_selected_endocytic_genes[, c(1:20)] %>%
  filter_if(is.numeric, any_vars(. <= -0.323))


setwd(file.path("~/Desktop/PJATK project"))
dfDEG.file.cases.all <- unname(unlist(read.table("dfDEG.file.cases.all.txt")))

Volcano_Plot <- function(dataset, FDR_cutoff = 0.05, UP = 0.50, DOWN = -0.50) {
  directory <- setwd(file.path("~/Desktop/PJATK project/dataDEG"))
  dfDEG <- read.table(paste("dfDEG_TCGA-", dataset, ".txt", sep = ""), header = TRUE)
  dfDEG$Threshold <- with(dfDEG, ifelse(dfDEG$logFC >= UP & dfDEG$FDR < FDR_cutoff, "Upregulated",
                                        ifelse(dfDEG$logFC <= DOWN & dfDEG$FDR < FDR_cutoff,
                                               "Downregulated", "Not significant"))) 
  
  ESCRT_list <- c("HGS", "STAM", "STAM2",
                  "TSG101", "VPS28", "VPS37A", "VPS37B", "VPS37C", "VPS37D", 
                  "UBAP1", "MVB12A", "MVB12B",
                  "VPS36", "SNF8", "VPS25",
                  "CHMP6", "CHMP4A", "CHMP4B", "CHMP4C",
                  "CHMP3", "CHMP2A", "CHMP2B", "CHMP5", "PDCD6IP",
                  "IST1", "CHMP1", "VPS4A", "VPS4B")
  ESCRT <- filter(dfDEG, Symbol %in% ESCRT_list)
  
  vp <- ggplot(data = dfDEG, 
               mapping = aes(x = logFC, y = -log10(FDR), colour = Threshold)) +
    scale_color_manual(values = c("dodgerblue", "gold", "deeppink2")) +
    geom_point(alpha = 0.4, size = 1.0) + 
    xlim(c(-3.5, 3.5)) + ylim(c(0, 30)) + labs(color = "Expression pattern") +
    geom_point(data=ESCRT, colour="black") +
    ggtitle(paste("All samples for ", dataset, " data set - Tumor vs. Normal", sep = "")) + 
    xlab("log2FoldChange") + ylab("-log10(FDR)") +
    theme(plot.title = element_text(hjust = 0.5))
  vp + geom_text_repel(data = ESCRT, aes(label = ESCRT$Symbol), colour = "black", size = 4)
}


for (df in dfDEG.file.cases.all) {
  i <- gsub(pattern = "dfDEG_TCGA-", replacement = "", x = df)
  print(Volcano_Plot(dataset = i))
  ggsave(Volcano_Plot(dataset = i), 
         file=paste0("plot_dfDEG_TCGA_", i,".png"), width = 14, height = 10, units = "cm")
}