rm(list = ls())

# this script requires R 3.4.4 because of the VennDiagram package

############
### LIBS ###


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
  library(VennDiagram)
})

###########
### FUN ###

draw.2.Venn <- function (area1, area2, intersection, 
                         stage1 = "Early stages", stage2 = "Advanced stages") {
  draw.pairwise.venn(area1 = length(area1),
                     area2 = length(area2),
                     cross.area = length(intersection),
                     # c(paste0(stage1), paste0(stage2)), 
                     fill = c("dodgerblue1", "deeppink2"),
                     fontfamily = "Arial", cat.fontfamily = "Arial",
                     cex = 0, # 4
                     # cex = 4,
                     # cat.cex = 4, #
                     # cat.pos = c(30, 215), #
                     # cat.dist = .05, #
                     scaled = TRUE)
}

############
### MAIN ###

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

rm(DGE_StageIandII, DGE_StageIIIandIV, f, temp.file.name, file.list)

flog.debug("Identify trafficking-related genes in early and late stages of cancer progression")

early <- filter(DGE_Earlystages, DGE_Earlystages$Symbol %in% endocytic_genes$`Gene Symbol`)
late <- filter(DGE_Latestages, DGE_Latestages$Symbol %in% endocytic_genes$`Gene Symbol`)

rm(endocytic_genes)

flog.debug("Genes downreg in early stages of cancer progression")
early_up = filter(early, early$Threshold == "Upregulated")
early_down = filter(early, early$Threshold == "Downregulated")
early = rbind(early_up, early_down)
saveRDS(early, file.path(proj.dir, spec.dir, "stage_early.RDS"))
saveRDS(early_down, file.path(proj.dir, spec.dir, "stage_early_down.RDS"))
saveRDS(early_up, file.path(proj.dir, spec.dir, "stage_early_up.RDS"))

flog.debug("Genes downreg in late stages of cancer progression")
late_up = filter(late, late$Threshold == "Upregulated")
late_down = filter(late, late$Threshold == "Downregulated")
late = rbind(late_up, late_down)
saveRDS(late, file.path(proj.dir, spec.dir, "stage_late.RDS"))
saveRDS(late_down, file.path(proj.dir, spec.dir, "stage_late_down.RDS"))
saveRDS(late_up, file.path(proj.dir, spec.dir, "stage_late_up.RDS"))

# rm(early, late)

flog.debug("Genes at the union and intersection between stages")
intersect <- intersect(early$Symbol, late$Symbol)
unique_4_early <- setdiff(early$Symbol, intersect)
cat("number of unique genes in early stages: ", length(unique_4_early))
unique_4_late <- setdiff(late$Symbol, intersect)
cat("number of unique genes in advanced stages: ", length(unique_4_late))

flog.debug("Venn diagram")
draw.2.Venn(late$Symbol, early$Symbol, intersect)


###################
### SessionInfo ###
###################

sessionInfo()