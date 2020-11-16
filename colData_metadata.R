# colData_metadata.R

library(tidyverse)
library(stringr)
library(hash)
library(GEOquery)


# get matadata info from NCBI/GEO -----------------------------------------

download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Run_Members.tab",
              destfile = "SRA_Run_Members.tab")
srr_srx <- read_tsv("SRA_Run_Members.tab")
gse_srp <- read_tsv("GSE_SRP.txt")  # written by hand
srp_gse_hash <- hash(gse_srp$SRP, gse_srp$GSE)

toBeDownload <- filter(srr_srx, Study %in% gse_srp$SRP)
# or use SRX in the pData from getGEO
unique.default(toBeDownload$Study)  # 7 GSEs, none are missed

# add gse info to the 'to be download' table
gse <- c()
for(srp in toBeDownload$Study){
  gse <-c(gse, srp_gse_hash[[srp]])
}
toBeDownload <- cbind(toBeDownload, gse)
saveRDS(toBeDownload, file = "to_be_downloaded_from_GEO.rds")
# to be used by 


### get the latest phenotype data from GEO

esets <- list()
for(gse in gse_srp$GSE){
  esets <- c(esets, getGEO(gse))
}
saveRDS(esets, file = "esets_from_geo.rds")


# importing metadata ---------------------------------------------------------

## GSEs
setwd("/binf-isilon/alab/students/vrw936/Master_Thesis/MyCode")
esets <- readRDS("esets_from_geo.rds")  # GSE122438 should be dumped

## check sample sizes across genders and diagnosis
# 1
View(pData(esets$GSE125583_series_matrix.txt.gz))
names(pData(esets$GSE125583_series_matrix.txt.gz))
pheno_data <- select(pData(esets$GSE125583_series_matrix.txt.gz), `case id:ch1`,
                     `age:ch1`, `characteristics_ch1.4`, `diagnosis:ch1`)
pheno_data <- unique.data.frame(pheno_data)
table(select(dplyr::filter(pheno_data, 
                           `age:ch1` >= 65), 
             characteristics_ch1.4, `diagnosis:ch1`))

# 2
View(pData(esets$GSE125050_series_matrix.txt.gz))
names(pData(esets$GSE125050_series_matrix.txt.gz))
pheno_data <- select(pData(esets$GSE125050_series_matrix.txt.gz), characteristics_ch1,
                     `expired_age:ch1`, `Sex:ch1`, `diagnosis:ch1`)
pheno_data <- unique.data.frame(pheno_data)
table(select(dplyr::filter(pheno_data, 
                           `expired_age:ch1` >= 65), 
             `Sex:ch1`, `diagnosis:ch1`))

# 3
View(pData(esets$GSE110731_series_matrix.txt.gz))
names(pData(esets$GSE110731_series_matrix.txt.gz))
pheno_data <- select(pData(esets$GSE110731_series_matrix.txt.gz),
                     `age:ch1`, `Sex:ch1`, `braak stage:ch1`)
# pheno_data <- unique.data.frame(pheno_data)
table(select(dplyr::filter(pheno_data, 
                           `age:ch1` >= 65), 
             `Sex:ch1`, `braak stage:ch1`))

# 4 GSE104704, no AD diagnosis or gender info shown
View(pData(esets$GSE104704_series_matrix.txt.gz))
names(pData(esets$GSE104704_series_matrix.txt.gz))
# pheno_data <- select(pData(esets$GSE104704_series_matrix.txt.gz),
#                      `age (years):ch1`, `Sex:ch1`, `braak stage:ch1`)
# pheno_data <- unique.data.frame(pheno_data)
# table(select(dplyr::filter(pheno_data, 
#                            `age:ch1` >= 65), 
#              `Sex:ch1`, `braak stage:ch1`))

# 5
View(pData(esets$GSE95587_series_matrix.txt.gz))
names(pData(esets$GSE95587_series_matrix.txt.gz))
pheno_data <- select(pData(esets$GSE95587_series_matrix.txt.gz), 
                     `age at death:ch1`, `Sex:ch1`, `diagnosis:ch1`)
# pheno_data <- unique.data.frame(pheno_data)
table(select(dplyr::filter(pheno_data, 
                           `age at death:ch1` >= 65), 
             `Sex:ch1`, `diagnosis:ch1`))

# 6 GSE53697 no gender or age info found
View(pData(esets$GSE53697_series_matrix.txt.gz))
names(pData(esets$GSE53697_series_matrix.txt.gz))
# pheno_data <- select(pData(esets$GSE53697_series_matrix.txt.gz), 
#                      `age at death:ch1`, `Sex:ch1`, `disease status:ch1`)
# pheno_data <- unique.data.frame(pheno_data)
# table(select(dplyr::filter(pheno_data, 
#                            `age at death:ch1` >= 65), 
#              `Sex:ch1`, `disease status:ch1`))

## MayoRNAseq
setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/MayoRNAseq/Gene Expression/Gene Expression (RNA seq)/cerebellum")
MayoRNAseq_cerebellum_cov <- read_csv("MayoRNAseq_RNAseq_CER_covariates.csv")
MayoRNAseq_cerebellum_qc <- read_tsv("MayoRNASeq_RNASeq_CBE_QCdetails.txt")
names(MayoRNAseq_cerebellum_cov)
MayoRNAseq_cerebellum_cov$AgeAtDeath <- as.numeric(gsub("(\\d{2})_or_above", 
                                                        "\\1", 
                                                        MayoRNAseq_cerebellum_cov$AgeAtDeath))
table(select(dplyr::filter(MayoRNAseq_cerebellum_cov,
                           AgeAtDeath >= 65),
             Sex, Diagnosis))
# AD, Control, Pathologic Aging, PSP? what do the last mean?
# obviouly only AD & Control are included in the number shown in the description
# samples that did not pass QC are excluded and described in the QC file

setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/MayoRNAseq/Gene Expression/Gene Expression (RNA seq)/temporal cortex")
MayoRNAseq_tcx_cov <- read_csv("MayoRNAseq_RNAseq_TCX_covariates.csv")
MayoRNAseq_tcx_qc <- read_tsv("MayoRNAseq_RNAseq_TCX_QCdetails.txt")
names(MayoRNAseq_tcx_cov)
MayoRNAseq_tcx_cov$AgeAtDeath <- as.numeric(gsub("(\\d{2})_or_above", 
                                                        "\\1", 
                                                 MayoRNAseq_tcx_cov$AgeAtDeath))
table(select(dplyr::filter(MayoRNAseq_tcx_cov,
                           AgeAtDeath >= 65),
             Gender, Diagnosis))
# AD, Control, Pathologic Aging, PSP? what do the last mean?
# obviouly only AD & Control are included in the number shown in the description
# samples that did not pass QC are excluded and described in the QC file


## MCADGS, no diagnosis info found??? see the paper
setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/MCADGS/MayoPilotRNAseq/Metadata")
MCADGS <- read_csv("MayoPilotRNAseq_RNAseq_TCX_AD_covariates.csv")
table(select(dplyr::filter(MCADGS,
                           `Age at death` >= 65),
             Sex))

## MSBB
setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/MSBB/Metadata")
MSBB_clinical <- read_csv("MSBB_clinical.csv")
MSBB_RNAseq_cov <- read_csv("MSBB_RNAseq_covariates_November2018Update.csv")
MSBB_clinical$diagnosis <- rep("", dim(MSBB_clinical)[1])
MSBB_clinical$AOD <- as.numeric(gsub("(\\d{2})\\+", 
                                                 "\\1", 
                                            MSBB_clinical$AOD))

MSBB_clinical$diagnosis[MSBB_clinical$bbscore <= 2.5] <- "Control"
MSBB_clinical$diagnosis[MSBB_clinical$bbscore >= 4.0] <- "AD"
MSBB_clinical$diagnosis[(MSBB_clinical$bbscore < 4.0) && (MSBB_clinical$bbscore > 2.5)] <- "out"
table(select(dplyr::filter(MSBB_clinical,
                           AOD >= 65),
             SEX, diagnosis))

BA <- unique.data.frame(select(MSBB_RNAseq_cov, individualIdentifier, BrodmannArea))
BA <- na.omit(BA)
age <- hash(MSBB_clinical$individualIdentifier, MSBB_clinical$AOD)
sex <- hash(MSBB_clinical$individualIdentifier, MSBB_clinical$SEX)
diagnosis <- hash(MSBB_clinical$individualIdentifier, MSBB_clinical$diagnosis)

BA$age <- rep(0, dim(BA)[1])
BA$sex <- rep("", dim(BA)[1])
BA$diagnosis <- rep("", dim(BA)[1])

for(i in 1:dim(BA)[1]){
  BA$age[i] <- age[[BA$individualIdentifier[i]]]
  BA$sex[i] <- sex[[BA$individualIdentifier[i]]]
  BA$diagnosis[i] <- diagnosis[[BA$individualIdentifier[i]]]
}

table(select(dplyr::filter(BA,
                           age >= 65),
             sex, diagnosis, BrodmannArea))

## ROSMAP
setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/ROSMAP/Metadata")
ROSMAP_RNAseq <- read_csv("ROSMAP_assay_RNAseq_metadata.csv")
ROSMAP_microRNAarray <- read_csv("ROSMAP_arraymiRNA_covariates.csv")
ROSMAP_biospecimen <- read_csv("ROSMAP_biospecimen_metadata.csv")
ROSMAP_clinical <- read_csv("ROSMAP_Clinical_2019-05_v3.csv")

specimen_rnaseq <- filter(ROSMAP_biospecimen, notes == "geneExpression (rnaSeq)")
specimen_rnaseq <- na.omit(select(specimen_rnaseq, individualID, tissue))
specimen_rnaseq <- unique.data.frame(specimen_rnaseq)
ROSMAP_clinical$diagnosis <- rep("", dim(ROSMAP_clinical)[1])
for(i in 1:dim(ROSMAP_clinical)[1]){
  if (is.na(ROSMAP_clinical$braaksc[i])) {
    ROSMAP_clinical$diagnosis[i] <- "out"
  }
  else if (ROSMAP_clinical$braaksc[i] <= 2.5) {
    ROSMAP_clinical$diagnosis[i] <- "Control"
  }
  else if (ROSMAP_clinical$braaksc[i] >= 4.0) {
    ROSMAP_clinical$diagnosis[i] <- "AD"
  }
  else{ROSMAP_clinical$diagnosis[i] <- "out"}
}
ROSMAP_clinical$age_at_visit_max <- as.numeric(gsub("(\\d{2})\\+", 
                                                    "\\1", 
                                                    ROSMAP_clinical$age_at_visit_max))
  
ROSMAP_clinical$msex <- gsub("1", "M", ROSMAP_clinical$msex)
ROSMAP_clinical$msex <- gsub("0", "F", ROSMAP_clinical$msex)

age <- hash(ROSMAP_clinical$individualID, ROSMAP_clinical$age_at_visit_max)
sex <- hash(ROSMAP_clinical$individualID, ROSMAP_clinical$msex)
diagnosis <- hash(ROSMAP_clinical$individualID, ROSMAP_clinical$diagnosis)

specimen_rnaseq$age <- rep(0, dim(specimen_rnaseq)[1])
specimen_rnaseq$sex <- rep("", dim(specimen_rnaseq)[1])
specimen_rnaseq$diagnosis <- rep("", dim(specimen_rnaseq)[1])

for(i in 1:dim(specimen_rnaseq)[1]){
  specimen_rnaseq$age[i] <- age[[specimen_rnaseq$individualID[i]]]
  specimen_rnaseq$sex[i] <- sex[[specimen_rnaseq$individualID[i]]]
  specimen_rnaseq$diagnosis[i] <- diagnosis[[specimen_rnaseq$individualID[i]]]
}

table(select(dplyr::filter(specimen_rnaseq,
                           age >= 65),
             diagnosis, tissue))




# make colData for each GSE to test DE analysis ---------------------------

## prepare the SAMN - SRR table
setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna")
GSEs<- grep(pattern = "^GSE\\d+", x = list.files(), value = TRUE)

srr_srx_samn <- readRDS("/binf-isilon/alab/students/vrw936/Master_Thesis/MyCode/to_be_downloaded_from_GEO.txt")
srrs <- list.files(file.path(".", GSEs, "trimmed", "salmon"))
biosample_run <- srr_srx_samn[srr_srx_samn$Run %in% srrs, c("BioSample", "Run")]

## functionalise
add_SAMN_SRR_GSE <- function(colData, biosample_run, gse){
  # a function adding SAMN, SRR, GSE 
  # to each line of a nearly ready-made colData table
  # and save the colData to an RDS file
  
  colData$GSE <- gse
  colData$BioSample <- gsub("^.+(SAMN[0-9]+)$", "\\1", colData$BioSample)
  colData <- inner_join(biosample_run, colData)
  row.names(colData) <- colData$Run
  colData <- as.matrix(colData)
  
  saveRDS(colData, paste0("colData_", gse, ".rds"))
  
  return(colData)
}


## make colData for each GSE

# GSE110731
colData_gse110731 <- pData(esets$GSE110731_series_matrix.txt.gz)
colData_gse110731 <- select(colData_gse110731, BioSample = relation, 
                            age = `age:ch1`, sex = `Sex:ch1`,
                            braak = `braak stage:ch1`, pmi = `pmi:ch1`)

colData_gse110731$condition = ""
colData_gse110731[colData_gse110731$braak %in% c(1, 2, "1,2"), "condition"] = "ctrl"
colData_gse110731[colData_gse110731$braak %in% c(4, 5, 6, "5,6"), "condition"] = "AD"
colData_gse110731$age <- as.numeric(colData_gse110731$age)
colData_gse110731 <- filter(colData_gse110731, (condition != "" & age >= 65))

colData_gse110731 <- add_SAMN_SRR_GSE(colData_gse110731, biosample_run, "GSE110731")


# GSE53697
colData_gse53697 <- pData(esets$GSE53697_series_matrix.txt.gz)
colData_gse53697 <- select(colData_gse53697, BioSample = relation, 
                           condition = `disease status:ch1`)
colData_gse53697$age <- ""
colData_gse53697$sex <- ""
colData_gse53697$braak <- ""
colData_gse53697$pmi <- ""

colData_gse53697[colData_gse53697$condition == "advanced Alzheimer's Disease", "condition"] = "AD"
colData_gse53697[colData_gse53697$condition == "control", "condition"] = "ctrl"

colData_gse53697 <- add_SAMN_SRR_GSE(colData_gse53697, biosample_run, "GSE53697")
# PROBLEM: after this inner_join, only 1 of the 9 AD samples are left
# but all SRR's salmon quant result exist. maybe the SAMN_SRR_GSE table is not up-to-date


# GSE104704
colData_gse104704 <- pData(esets$GSE104704_series_matrix.txt.gz)
colData_gse104704 <- select(colData_gse104704, BioSample = relation, 
                            age = `age (years):ch1`, condition = `study group:ch1`)

colData_gse104704$age <- as.numeric(colData_gse104704$age)
colData_gse104704 <- filter(colData_gse104704, age >= 65)
colData_gse104704$condition <- sub("Aged, diseased", "AD", colData_gse104704$condition)
colData_gse104704$condition <- sub("Old", "ctrl", colData_gse104704$condition)

colData_gse104704$braak <- ""
colData_gse104704$sex <- ""
colData_gse104704$pmi <- ""

colData_gse104704 <- add_SAMN_SRR_GSE(colData_gse104704, biosample_run, "GSE104704")

# GSE125583
colData_gse125583 <- pData(esets$GSE125583_series_matrix.txt.gz)
colData_gse125583 <- select(colData_gse125583, BioSample = relation, 
                            age = `age:ch1`, sex = `Sex:ch1`,
                            braak = `braak.score:ch1`, 
                            condition = `diagnosis:ch1`)

colData_gse125583 <- filter(colData_gse125583, !str_detect(BioSample, "Alternative"))
# remove samples that are repeated in GSE95587

colData_gse125583$age <- as.numeric(colData_gse125583$age) # without this, age above 100 will be lost
colData_gse125583 <- filter(colData_gse125583, age >= 65)

colData_gse125583$pmi <- ""
colData_gse125583$condition <- sub("Alzheimer's disease", "AD", colData_gse125583$condition)
colData_gse125583$condition <- sub("control", "ctrl", colData_gse125583$condition)

colData_gse125583 <- add_SAMN_SRR_GSE(colData_gse125583, biosample_run, "GSE125583")
# PROBLEM: some ctrl has braak III or IV, some AD has braak I/II/III

# GSE125050
colData_gse125050 <- pData(esets$GSE125050_series_matrix.txt.gz)
colData_gse125050 <- select(colData_gse125050, BioSample = relation, 
                            age = `expired_age:ch1`, sex = `Sex:ch1`,
                            braak = `braak.score:ch1`, pmi = `pmi:ch1`,
                            condition = `diagnosis:ch1`, `cell type:ch1`,
                            patient = `inventory_patient_id:ch1`)

# PROBLEM: originally set braak I,II,III,IV as control, V,VI as AD
colData_gse125050$age <- as.numeric(colData_gse125050$age)
colData_gse125050 <- filter(colData_gse125050, age >= 65)
colData_gse125050$condition <- sub("Control", "ctrl", colData_gse125050$condition)
# PROBLEM: gender of some samples missing, need to be predicted.
colData_gse125050 <- add_SAMN_SRR_GSE(colData_gse125050, biosample_run, "GSE125050")



# GSE95587
colData_gse95587 <- pData(esets$GSE95587_series_matrix.txt.gz)

colData_gse95587$relation <- as.character(colData_gse95587$relation)
bs_tmp <- as.character(colData_gse95587[str_detect(colData_gse95587$relation, "Alternative"), "relation.1"])
colData_gse95587[str_detect(colData_gse95587$relation, "Alternative"), "relation"] <- bs_tmp
# for sample sequenced in both GSE95587 and GSE125583, keep the data in GSE95587

colData_gse95587 <- select(colData_gse95587, BioSample = relation, 
                           age = `age at death:ch1`, sex = `Sex:ch1`,
                           braak = `braak score:ch1`, pmi = `post_mortem_interval:ch1`,
                           condition = `diagnosis:ch1`, rin = `rin_score:ch1`)


colData_gse95587$condition <- sub("Alzheimer's disease", "AD", colData_gse95587$condition)
colData_gse95587$condition <- sub("control", "ctrl", colData_gse95587$condition)

colData_gse95587$age <- as.numeric(colData_gse95587$age)  # this 2 lines can be functionalise
colData_gse95587 <- filter(colData_gse95587, age >= 65)

colData_gse95587 <- add_SAMN_SRR_GSE(colData_gse95587, biosample_run, "GSE95587")

## make colData for each AMP-AD datasets

# ROSMAP (2020.5.14: actually only include FL fastq files;)
setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/ROSMAP/Metadata")
ROSMAP_RNAseq <- read_csv("ROSMAP_assay_RNAseq_metadata.csv")
ROSMAP_biospecimen <- read_csv("ROSMAP_biospecimen_metadata.csv")
ROSMAP_clinical <- read_csv("ROSMAP_Clinical_2019-05_v3.csv")

sample_id <- list.files("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/ROSMAP/fastq/")
sample_id <- sample_id[str_detect(sample_id, ".fastq.gz")]
sample_id <- unique.default(gsub("(^.+)_[1|2]\\.fastq\\.gz", "\\1", sample_id))

ROSMAP_RNAseq <- filter(ROSMAP_RNAseq, specimenID %in% sample_id)
ROSMAP_RNAseq <- select(ROSMAP_RNAseq, specimenID, RIN, libraryBatch, sequencingBatch)
ROSMAP_RNAseq <- filter(ROSMAP_RNAseq, sequencingBatch != "nan")

ROSMAP_biospecimen <- filter(ROSMAP_biospecimen, notes == "geneExpression (rnaSeq)")
ROSMAP_biospecimen <- na.omit(select(ROSMAP_biospecimen, individualID, specimenID, tissue))

# prepare clinical (about criteria, only braak.) and inner_join the 3 table
# fastq files of 2 brain regions are missing, need to be requested 
# or generated from BAM by myself

ROSMAP_clinical$diagnosis <- rep("", dim(ROSMAP_clinical)[1])
for(i in 1:dim(ROSMAP_clinical)[1]){
  if (is.na(ROSMAP_clinical$braaksc[i])) {
    ROSMAP_clinical$diagnosis[i] <- "out"
  }
  else if (ROSMAP_clinical$braaksc[i] <= 2.5) {
    ROSMAP_clinical$diagnosis[i] <- "ctrl"
  }
  else if (ROSMAP_clinical$braaksc[i] >= 4.0) {
    ROSMAP_clinical$diagnosis[i] <- "AD"
  }
  else{ROSMAP_clinical$diagnosis[i] <- "out"}
}
ROSMAP_clinical$age_at_visit_max <- round(as.numeric(gsub("(\\d{2})\\+", 
                                                          "\\1", 
                                                          ROSMAP_clinical$age_at_visit_max)))

ROSMAP_clinical$msex <- gsub("1", "M", ROSMAP_clinical$msex)
ROSMAP_clinical$msex <- gsub("0", "F", ROSMAP_clinical$msex)

colData_ROSMAP <- inner_join(ROSMAP_RNAseq, ROSMAP_biospecimen, by = "specimenID")
colData_ROSMAP <- inner_join(colData_ROSMAP, ROSMAP_clinical, by = "individualID")

# check the data and remove redundancy
length(unique.default(colData_ROSMAP$specimenID))
length(colData_ROSMAP$specimenID)
tmp <- (unique.default(colData_ROSMAP$specimenID) == colData_ROSMAP$specimenID)
colData_ROSMAP$specimenID[c(321, 322, 323)]
colData_ROSMAP <- colData_ROSMAP[c(-321, -322),]

colData_ROSMAP <- select(colData_ROSMAP, BioSample = specimenID, Run = specimenID, 
                         age = age_at_visit_max, sex = msex, condition = diagnosis,
                         braak = braaksc, pmi, individualID, ceradsc, cogdx, dcfdx_lv,
                         educ, race, spanish, apoe_genotype, tissue, GSE = Study)
# make sure the colData can be merged with the colData of GSEs datasets

colData_ROSMAP <- filter(colData_ROSMAP, condition != "out", age >= 65)
row.names(colData_ROSMAP) <- colData_ROSMAP$Run
colData_ROSMAP <- as.matrix(colData_ROSMAP)

saveRDS(colData_ROSMAP, "/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/colData_ROSMAP.rds")


# MayoRNAseq

setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/MayoRNAseq/Gene Expression/Gene Expression (RNA seq)/cerebellum")
MayoRNAseq_cerebellum_cov <- read_csv("MayoRNAseq_RNAseq_CER_covariates.csv")

MayoRNAseq_cerebellum_cov$AgeAtDeath <- as.numeric(gsub("(\\d{2})_or_above", 
                                                        "\\1", 
                                                        MayoRNAseq_cerebellum_cov$AgeAtDeath))
colData_MayoRNAseq_CBE <- select(MayoRNAseq_cerebellum_cov, BioSample = SampleID, 
                                 Run = SampleID, age = AgeAtDeath, sex = Sex,
                                 condition = Diagnosis, braak = Braak, rin = RIN,
                                 pmi = PMI, apoe_genotype = ApoE, Tissue, Flowcell, 
                                 Source, Thal)
colData_MayoRNAseq_CBE$GSE = "MayoRNAseq"
colData_MayoRNAseq_CBE$condition <- sub("Control", "ctrl", colData_MayoRNAseq_CBE$condition)


colData_MayoRNAseq_CBE <- filter(colData_MayoRNAseq_CBE, age >= 65, condition == "AD" | condition == "ctrl")
# 278 -> 150
# AD, Control, Pathologic Aging, PSP? what does the last one mean?
# obviouly only AD & Control are included in the number shown in the description
# samples that did not pass QC are excluded and described in the QC file

rn_tmp <- colData_MayoRNAseq_CBE$Run
colData_MayoRNAseq_CBE <- as.matrix(colData_MayoRNAseq_CBE)
row.names(colData_MayoRNAseq_CBE) <- rn_tmp
saveRDS(colData_MayoRNAseq_CBE, "/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/colData_MayoRNAseq_CBE.rds")



setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/MayoRNAseq/Gene Expression/Gene Expression (RNA seq)/temporal cortex")
MayoRNAseq_tcx_cov <- read_csv("MayoRNAseq_RNAseq_TCX_covariates.csv")
MayoRNAseq_tcx_cov$AgeAtDeath <- as.numeric(gsub("(\\d{2})_or_above", 
                                                 "\\1", 
                                                 MayoRNAseq_tcx_cov$AgeAtDeath))

colData_MayoRNAseq_TCX <- select(MayoRNAseq_tcx_cov, BioSample = SampleID, 
                                 Run = SampleID, age = AgeAtDeath, sex = Gender,
                                 condition = Diagnosis, braak = Braak, rin = RIN,
                                 pmi = PMI, apoe_genotype = ApoE, Tissue, FLOWCELL, 
                                 Source, Thal)
colData_MayoRNAseq_TCX$GSE = "MayoRNAseq"
colData_MayoRNAseq_TCX$condition <- sub("Control", "ctrl", colData_MayoRNAseq_TCX$condition)


colData_MayoRNAseq_TCX <- filter(colData_MayoRNAseq_TCX, age >= 65, condition == "AD" | condition == "ctrl")
# 278 -> 152
# AD, Control, Pathologic Aging, PSP? what does the last one mean?
# obviouly only AD & Control are included in the number shown in the description
# samples that did not pass QC are excluded and described in the QC file

rn_tmp <- colData_MayoRNAseq_TCX$Run
colData_MayoRNAseq_TCX <- as.matrix(colData_MayoRNAseq_TCX)
row.names(colData_MayoRNAseq_TCX) <- rn_tmp
saveRDS(colData_MayoRNAseq_TCX, "/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/colData_MayoRNAseq_TCX.rds")


# MSBB
setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/MSBB/Metadata")
MSBB_clinical <- read_csv("MSBB_clinical.csv")
MSBB_RNAseq_cov <- read_csv("MSBB_RNAseq_covariates_November2018Update.csv")
MSBB_clinical$diagnosis <- rep("", dim(MSBB_clinical)[1])
MSBB_clinical$AOD <- as.numeric(gsub("(\\d{2})\\+", 
                                     "\\1", 
                                     MSBB_clinical$AOD))

MSBB_clinical$diagnosis[MSBB_clinical$bbscore <= 2.5] <- "ctrl"
MSBB_clinical$diagnosis[MSBB_clinical$bbscore >= 4.0] <- "AD"

MSBB_clinical <- filter(MSBB_clinical, AOD >= 65, diagnosis == "AD" | diagnosis == "ctrl")


MSBB_RNAseq_cov <- unique.data.frame(select(MSBB_RNAseq_cov, sampleIdentifier, 
                                            individualIdentifier, BrodmannArea,
                                            batch, RIN))

colData_MSBB <- inner_join(MSBB_RNAseq_cov, MSBB_clinical, by = "individualIdentifier")

length(unique.default(colData_MSBB$sampleIdentifier))
length(colData_MSBB$sampleIdentifier)
tmp <- (colData_MSBB$sampleIdentifier == unique.default(colData_MSBB$sampleIdentifier))

colData_MSBB <- select(colData_MSBB, BioSample = sampleIdentifier, Run = sampleIdentifier,
                       age = AOD, sex = SEX, condition = diagnosis, braak = bbscore,
                       rin = RIN, pmi = PMI, RACE, CDR, `NP.1`, PlaqueMean, batch,
                       BrodmannArea)
colData_MSBB$GSE <- "MSBB"

rn_tmp <- colData_MSBB$Run
colData_MSBB <- as.matrix(colData_MSBB)
row.names(colData_MSBB) <- rn_tmp

saveRDS(colData_MSBB, "/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/colData_MSBB.rds")

colData_MSBB_TL <- colData_MSBB[colData_MSBB[, "BrodmannArea"] %in% c("BM22","BM36"), ]
colData_MSBB_FL <- colData_MSBB[colData_MSBB[, "BrodmannArea"] %in% c("BM10","BM44"), ]

saveRDS(colData_MSBB_TL, "/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/colData_MSBB_TL.rds")
saveRDS(colData_MSBB_FL, "/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/colData_MSBB_FL.rds")




