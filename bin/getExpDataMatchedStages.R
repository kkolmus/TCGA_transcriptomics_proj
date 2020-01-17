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


flog.debug("Load information required for the analysis")

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl = TRUE)]
projects

CancerProject <- c("TCGA-COAD", "TCGA-READ")
DataDirectory <- "GDCdata"
platform <- "Illumina HiSeq"
FileNameData <- paste0(DataDirectory, platform,".rda")
file.type <- "results"

dataClin <- readRDS(file.path(proj.dir, "data", "ClinData.RDS"))
dataClin$ajcc_pathologic_stage <- as.factor(dataClin$ajcc_pathologic_stage)
dataClin$ajcc_pathologic_t <- as.factor(dataClin$ajcc_pathologic_t)

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
rm(samplesDown)


flog.debug("Download, normalize, filter and perform differential gene expression analysis")

DEGanalysis <- function(stages, title,
                        UP = 0.6, DOWN = 0.6, FDR_cutoff = 0.05,
                        PreProc_cor.cut = 0.6, 
                        Norm_method = "gcContent", 
                        Filt_method = "quantile", 
                        Filt_qnt.cut = 0.25,
                        DEA_batch.factor = "Plate", 
                        DEA_method = "glmLRT") {
  
  title = str_replace_all(title, fixed(" "), "")
  
  # Select only samples with clinical information
  flog.debug("Select only samples with clinical information")
  for(s in stages) {
    Stage_NT <- filter(dataClin, dataClin$ajcc_pathologic_stage == s)
    SampleNT_Stage <- filter(Stage_NT, Stage_NT$bcr_patient_barcode %in% Sample.Short)
    SampleNT_Stage <- pmatch(SampleNT_Stage$bcr_patient_barcode, SampleNT)
    SampleNT_Stage <- SampleNT[SampleNT_Stage]
    
    SampleNT_Stage_final <- c(SampleNT_Stage_final, SampleNT_Stage)
    
    Stage_TP <- filter(dataClin, dataClin$ajcc_pathologic_stage == s)
    SampleTP_Stage <- filter(Stage_TP, Stage_TP$bcr_patient_barcode %in% Sample.Short)
    SampleTP_Stage <- pmatch(SampleTP_Stage$bcr_patient_barcode, SampleTP)
    SampleTP_Stage <- SampleTP[SampleTP_Stage]
    
    SampleTP_Stage_final <- c(SampleTP_Stage_final, SampleTP_Stage)
  }
  
  # Query platform Illumina HiSeq to download samples
  flog.debug("Query GDC")
  assign("queryDown", 
         GDCquery(project = CancerProject,
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = platform,
                  file.type = file.type,
                  barcode = c(SampleTP_Stage_final, SampleNT_Stage_final),
                  experimental.strategy = "RNA-Seq",
                  legacy = TRUE))
  saveRDS(queryDown, file.path(proj.dir, "data_matchedStages", paste0("GDCquery_", title, ".RDS")))
  # Download samples
  flog.debug("Download samples")
  tryCatch(GDCdownload(query = queryDown,
                       method = "api", 
                       files.per.chunk = 20,
                       directory = file.path(proj.dir, "data_matchedStages", "GDCdata")),
           error = function(e) GDCdownload(query = queryDown,
                                           method = "client", 
                                           files.per.chunk = 20,
                                           directory = file.path(proj.dir, "data_matchedStages", "GDCdata")))
  # Prepare samples for analysis
  flog.debug("Prepare GDC data")
  FileNameData <- paste0("GDCdata", str_replace_all(platform, fixed(" "), ""),".rda")
  dataPrep <- GDCprepare(query = queryDown, save = TRUE,
                         directory = file.path(proj.dir, "data_matchedStages", "GDCdata"), 
                         save.filename = FileNameData)
  saveRDS(dataPrep, file.path(proj.dir, "data_matchedStages", paste0("GDCprepare_", title, ".RDS")))
  # Samples preprocessing
  flog.debug("Perform intensity correlation")
  dataPreProc <- TCGAanalyze_Preprocessing(object = dataPrep, cor.cut = PreProc_cor.cut)
  saveRDS(dataPreProc, file.path(proj.dir, "data_matchedStages", paste0("dataPreProc_", title, ".RDS")))
  # Samples normalisation
  flog.debug("Perform normalization")
  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPreProc,
                                        geneInfo = geneInfo,
                                        method = Norm_method)
  saveRDS(dataNorm, file.path(proj.dir, "data_matchedStages", paste0("dataNorm_", title, ".RDS")))
  # Data filtering
  flog.debug("Perform data filtering based on threshold defined quantile mean across all samples")
  dataFilt <- assign(paste0("dataFilt_", title), 
                     TCGAanalyze_Filtering(tabDF = dataNorm, 
                                           method = Filt_method, 
                                           qnt.cut =  Filt_qnt.cut))
  dataFilt_transposed <- t(dataFilt)
  dataFilt_transposed <- as.data.frame(dataFilt_transposed)
  dataFilt_transposed <- rownames_to_column(dataFilt_transposed, "barcode_id")
  # Selecting filtered data based on sample type
  F_SampleNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = "NT")
  F_SampleTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = "TP")
  # normal tissue
  dataFilt_transposed_NT <- filter(dataFilt_transposed, 
                                   dataFilt_transposed$barcode_id %in% F_SampleNT)
  dataFilt_transposed_NT <- column_to_rownames(dataFilt_transposed_NT, "barcode_id")
  dataFilt_transposed_NT <- t(dataFilt_transposed_NT)
  dataFilt_transposed_NT <- as.data.frame(dataFilt_transposed_NT)
  dataFilt_transposed_NT <- rownames_to_column(dataFilt_transposed_NT, "Symbol")
  # tumour tissue
  dataFilt_transposed_TP <- filter(dataFilt_transposed, 
                                   dataFilt_transposed$barcode_id %in% F_SampleTP)
  dataFilt_transposed_TP <- column_to_rownames(dataFilt_transposed_TP, "barcode_id")
  dataFilt_transposed_TP <- t(dataFilt_transposed_TP)
  dataFilt_transposed_TP <- as.data.frame(dataFilt_transposed_TP)
  dataFilt_transposed_TP <- rownames_to_column(dataFilt_transposed_TP, "Symbol")
  # combine dataframes
  flog.debug("Prepare dataframe with filtered data for patients included in the analysis")
  patients <<- assign(paste0("dataFilt_patients_", title), 
                      merge(dataFilt_transposed_NT, 
                            dataFilt_transposed_TP, by = "Symbol"))
  saveRDS(patients,
          file.path(proj.dir, "data_matchedStages", paste0("patients_", title,".RDS")))
  # Differential gene expression analysis 
  flog.debug("Perform differential gene expression analysis")
  dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,F_SampleNT], mat2 = dataFilt[,F_SampleTP],
                              Cond1type = "Normal", Cond2type = "Tumor",
                              batch.factors = DEA_batch.factor,
                              fdr.cut = 1, logFC.cut = 0,
                              method = DEA_method)
  
  dataDEGs <- tibble::rownames_to_column(dataDEGs, var = "Symbol")
  dataDEGs <- mutate_at(dataDEGs, vars(-Symbol), funs(as.numeric(.)))
  dataDEGs$Threshold <- with(dataDEGs, ifelse(dataDEGs$logFC >= UP & 
                                                dataDEGs$FDR < FDR_cutoff, 
                                              "Upregulated",
                                              ifelse(dataDEGs$logFC <= DOWN & 
                                                       dataDEGs$FDR < FDR_cutoff, 
                                                     "Downregulated", "Not significant")))
  dataDEGs <<- assign(paste0("dataDEGs_", title), dataDEGs)
  saveRDS(dataDEGs,
          file.path(proj.dir, "data_matchedStages", paste0("dataDEG_", title,".RDS")))
}


flog.debug("Perform analysis for the early stages of carcinogenesis")

SampleNT_Stage_final <- c()
SampleTP_Stage_final <- c()

DEGanalysis(stages = c("Stage I", "Stage II"), title = "Stage I and II")

SampleNT_Stage_final <- c()
SampleTP_Stage_final <- c()

DEGanalysis(stages = c("Stage I", "Stage IA", 
                       "Stage II", "Stage IIA", "Stage IIB", "Stage IIC"), 
            title = "Early stages")


flog.debug("Perform analysis for the late stages of carcinogenesis")

SampleNT_Stage_final <- c()
SampleTP_Stage_final <- c()

DEGanalysis(stages = c("Stage III", "Stage IV"), title = "Stages III and IV")

SampleNT_Stage_final <- c()
SampleTP_Stage_final <- c()

DEGanalysis(stages = c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", 
                       "Stage IV", "Stage IVA", "Stage IVB"), 
            title = "Late stages")