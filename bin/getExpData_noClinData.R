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

data.dir <- file.path(proj.dir, "data_matched_noClinInfo")
setwd(data.dir)

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

# MatchedCoupledSampleTypes <- TCGAquery_MatchedCoupledSampleTypes(samplesDown, c("NT","TP"))
# SampleNT <- TCGAquery_SampleTypes(barcode = MatchedCoupledSampleTypes, typesample = "NT")
# SampleTP <- TCGAquery_SampleTypes(barcode = MatchedCoupledSampleTypes, typesample = "TP")

SampleNT <- TCGAquery_SampleTypes(barcode = samplesDown, typesample = "NT")
SampleTP <- TCGAquery_SampleTypes(barcode = samplesDown, typesample = "TP")

rm(samplesDown)

flog.debug("Download, normalize, filter and perform differential gene expression analysis")

DEGanalysis <- function(SampleNT, SampleTP,
                        UP = 0.6, DOWN = 0.6, FDR_cutoff = 0.05,
                        PreProc_cor.cut = 0.6, 
                        Norm_method = "gcContent", # "geneLength"
                        Filt_method = "quantile", 
                        Filt_qnt.cut = 0.25,
                        DEA_batch.factor = "Plate", 
                        DEA_method = "glmLRT") {
  
  # Query platform Illumina HiSeq to download samples
  flog.debug("Query GDC")
  assign("queryDown", 
         GDCquery(project = CancerProject,
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = platform,
                  file.type = file.type,
                  barcode = c(SampleNT, SampleTP),
                  experimental.strategy = "RNA-Seq",
                  legacy = TRUE))
  saveRDS(queryDown, file.path(data.dir, "GDCquery.RDS"))
  # Download samples
  flog.debug("Download samples")
  tryCatch(GDCdownload(query = queryDown,
                       method = "api", 
                       files.per.chunk = 20,
                       directory = file.path(data.dir, "GDCdata")),
           error = function(e) GDCdownload(query = queryDown,
                                           method = "client", 
                                           files.per.chunk = 20,
                                           directory = file.path(data.dir, "GDCdata")))
  # Prepare samples for analysis
  flog.debug("Prepare GDC data")
  FileNameData <- paste0("GDCdata", str_replace_all(platform, fixed(" "), ""),".rda")
  dataPrep <- GDCprepare(query = queryDown, save = TRUE,
                         directory = file.path(data.dir, "GDCdata"), 
                         save.filename = FileNameData)
  saveRDS(dataPrep, file.path(data.dir, "GDCprepare.RDS"))
  # Samples preprocessing
  flog.debug("Perform intensity correlation")
  dataPreProc <- TCGAanalyze_Preprocessing(object = dataPrep, cor.cut = PreProc_cor.cut)
  saveRDS(dataPreProc, file.path(data.dir, paste0("dataPreProc.RDS")))
  # Samples normalisation
  flog.debug("Perform normalization")
  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPreProc,
                                        geneInfo = geneInfo,
                                        method = Norm_method)
  saveRDS(dataNorm, file.path(data.dir, paste0("dataNorm.RDS")))
  # Data filtering
  flog.debug("Perform data filtering based on threshold defined quantile mean across all samples")
  dataFilt <- assign("dataFilt", 
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
  patients <<- assign("dataFilt_patients", 
                      merge(dataFilt_transposed_NT, 
                            dataFilt_transposed_TP, by = "Symbol"))
  saveRDS(patients,
          file.path(data.dir, paste0("patients.RDS")))
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
  dataDEGs <<- assign("dataDEGs", dataDEGs)
  saveRDS(dataDEGs, file.path(data.dir, "dataDEG.RDS"))
}


flog.debug("Perform analysis for the early stages of carcinogenesis")

DEGanalysis(SampleNT, SampleTP)
