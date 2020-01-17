rm(list = ls())

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("TCGAbiolinks")

suppressPackageStartupMessages({
  library(futile.logger)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(tibble)
  library(stringr)
})


flog.threshold(DEBUG)


flog.debug("Set project directory")

proj.dir = "~/Desktop/HT projects/TCGA_transcriptomics_proj"
ifelse(test = !dir.exists(file.path(proj.dir)), 
       yes = c(dir.create(file.path(proj.dir)), 
               setwd(file.path(proj.dir))), 
       no = "Folder exists")
getwd()


flog.debug("Set project directory and load data")

dataClin <- readRDS(file.path(proj.dir, "data_matchedSamples", "ClinData_matchedSamples.RDS"))
dataClin$ajcc_pathologic_stage <- as.factor(dataClin$ajcc_pathologic_stage)
levels <- dataClin$ajcc_pathologic_stage
dataClin$days_to_death <- as.numeric(dataClin$days_to_death)
early_stage <- filter(dataClin, 
                      dataClin$ajcc_pathologic_stage %in% 
                        c("Stage I", "Stage II", "Stage IIA"))
late_stage <- filter(dataClin, 
                     dataClin$ajcc_pathologic_stage %in% 
                       c("Stage III", "Stage IIIB", "Stage IIIC", "Stage IV", "Stage IVA"))

flog.debug("Load information required for the analysis")

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl = TRUE)]
projects

CancerProject <- c("TCGA-COAD", "TCGA-READ")
DataDirectory <- "GDCdata"
platform <- "Illumina HiSeq"
FileNameData <- paste0(DataDirectory, platform,".rda")
file.type <- "results"


flog.debug("Query GDC")

query <- GDCquery(project = CancerProject, 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = platform,
                  file.type = file.type,
                  experimental.strategy = "RNA-Seq",
                  legacy = TRUE)

samplesDown <- query$results[[1]]$cases

MatchedCoupledSampleTypes <- TCGAquery_MatchedCoupledSampleTypes(samplesDown, c("NT","TP"))
SampleNT <- TCGAquery_SampleTypes(barcode = MatchedCoupledSampleTypes, typesample = "NT")
SampleTP <- TCGAquery_SampleTypes(barcode = MatchedCoupledSampleTypes, typesample = "TP")
Sample.Short <- unique(substr(x = MatchedCoupledSampleTypes, start = 1, stop = 12))


flog.debug("Perform differential gene expression analysis and check patient survival")

# Inputs
PreProc_cor.cut = 0.6 # spearman correlation in samples by samples
Norm_method = "gcContent" # geneLength
Filt_method = "quantile" #  varFilter, filter1, filter2
Filt_qnt.cut = 0.25 # qnt.cut is threshold selected as mean for filtering
# Filt.eta= 0.05 # eta is a paramter for filter1. default eta = 0.05.
# Filt.foldChange	= 1 # foldChange is a paramter for filter2. default foldChange = 1.
DEA_batch.factor = "Plate"
DEA_method = "glmLRT"

dataFilt <- NULL
dataDEGs <- NULL
VPS37A <- NULL
VPS37B <- NULL


for (s in Sample.Short) {
  print(s)
  SampleNT_index = pmatch(s, SampleNT)
  SampleNT_pair = SampleNT[SampleNT_index]
  SampleTP_index = pmatch(s, SampleTP)
  SampleTP_pair = SampleTP[SampleTP_index]
  
  flog.debug("Query GDC")
  GDCquery_temp <- GDCquery(project = CancerProject,
                            data.category = "Gene expression",
                            data.type = "Gene expression quantification",
                            platform = platform,
                            file.type = file.type,
                            barcode = c(SampleNT_pair, SampleTP_pair),
                            experimental.strategy = "RNA-Seq",
                            legacy = TRUE)
  saveRDS(GDCquery_temp, 
          file.path(proj.dir, "data_matchedSamples", "GDCquery", paste0("GDCquery_", s, ".RDS")))
  
  flog.debug("Download samples")
  tryCatch(GDCdownload(query = GDCquery_temp,
                       method = "api", 
                       files.per.chunk = 20,
                       directory = file.path(proj.dir, "data", "GDCdata")),
           error = function(e) GDCdownload(query = GDCquery_temp,
                                           method = "client", 
                                           files.per.chunk = 20,
                                           directory = file.path(proj.dir, "data", "GDCdata")))
  
  flog.debug("Prepare GDC data")
  dataPrep <- GDCprepare(query = GDCquery_temp, save = TRUE,
                              directory = file.path(proj.dir, "data", "GDCdata"))
  saveRDS(dataPrep, 
          file.path(proj.dir, "data_matchedSamples", paste0("GDCprepare_", s, ".RDS")))
  
  flog.debug("Perform correlation")
  dataPreProc <- TCGAanalyze_Preprocessing(object = dataPrep, cor.cut = PreProc_cor.cut)
  saveRDS(dataPreProc, 
          file.path(proj.dir, "data_matchedSamples", paste0("dataPreProc_", s, ".RDS")))
  
  flog.debug("Perform normalization")
  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPreProc,
                                        geneInfo = geneInfo,
                                        method = Norm_method)
  saveRDS(dataNorm, file.path(proj.dir, "data_matchedSamples", paste0("dataNorm_", s, ".RDS")))
  
  flog.debug("Perform filtering")
  dataFilt[[s]] <- TCGAanalyze_Filtering(tabDF = dataNorm, 
                                         method = Filt_method, 
                                         qnt.cut =  Filt_qnt.cut)
  saveRDS(dataFilt[[s]], file.path(proj.dir, "data_matchedSamples", paste0("dataFilt_", s, ".RDS")))
  
  mat_temp <- as.data.frame(dataFilt[[s]])
  mat_temp <- rownames_to_column(mat_temp, var = "Symbol")
  
  flog.debug("Perform paired comparison of gene counts")
  # this cannot be done for a single patient
  # dataDEGs[[s]] <- TCGAanalyze_DEA(MAT = mat_temp,
  #                                  batch.factors = DEA_batch.factor,
  #                                  fdr.cut = 1, logFC.cut = 0,
  #                                  paired = TRUE,
  #                                  Condtypes = c(SampleNT_pair, SampleTP_pair),
  #                                  method = DEA_method)
  # 
  # dataDEGs[[s]] <- TCGAanalyze_DEA(MAT = mat_temp,
  #                                  batch.factors = DEA_batch.factor,
  #                                  Condtypes = c(SampleNT_pair, SampleTP_pair),
  #                                  fdr.cut = 1, logFC.cut = 0,
  #                                  method = DEA_method)
  
  normal_tissue = filter(mat_temp, Symbol == "VPS37A")[[2]]
  tumor_tissue = filter(mat_temp, Symbol == "VPS37A")[[3]]
  VPS37A[[s]] <- tumor_tissue/normal_tissue
  
  normal_tissue = filter(mat_temp, Symbol == "VPS37B")[[2]]
  tumor_tissue = filter(mat_temp, Symbol == "VPS37B")[[3]]
  VPS37B[[s]] <- tumor_tissue/normal_tissue
  
}

flog.debug("Combine expression and clinical data")

exp.fold.change <- as.data.frame(cbind(VPS37A, VPS37B))
exp.fold.change <- rownames_to_column(exp.fold.change, var = "bcr_patient_barcode")
dataClin <- full_join(dataClin, exp.fold.change, by = "bcr_patient_barcode")
# a single patient is missing data on Stage
dataClin <- filter(dataClin, dataClin$ajcc_pathologic_stage != "NA")

saveRDS(dataClin, file.path(proj.dir, "data_matchedSamples", paste0("dataClinExpData.RDS")))


library(survival)
library(survminer)

hist(dataClin$ajcc_pathologic_stage)
hist(dataClin$ajcc_pathologic_stage)

dataClin <- saveRDS(file.path(proj.dir, "data_matchedSamples", paste0("dataClinExpData.RDS")))

dataClin_subset <- dataClin[, c(2,3,8,9,10,11,20,24,25)]


