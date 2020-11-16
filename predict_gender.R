# predict the gender for each sample
library(tidyverse)
library(stringr)
library(DESeq2)



# get transcripts and genes on chromosome X and Y -------------------------

### load the annotation files

setwd("/binf-isilon/alab/students/vrw936/scratch/reference_ome_and_anno/")
gencode.v33.annotation <- read_tsv("gencode.v33.annotation.gtf",
                                   col_names = FALSE, skip = 5)

gencode.v33.metadata.HGNC <- read_tsv("gencode.v33.metadata.HGNC",
                                      col_names = FALSE)
# gencode.v33.metadata.EntrezGene <- read_tsv("gencode.v33.metadata.EntrezGene",
#                                             col_names = FALSE)

tx2gene <- select(gencode.v33.metadata.HGNC,
                  TXNAME = X1,
                  GENEID = X2)

### function


ChromAnno <- function(anno_gtf = gencode.v33.annotation, chroms = c("chrX", "chrY"), 
                      tx2gene_local = tx2gene){  # using tx2gene = tx2gene causes error
  # a function to extract transcripts on given chromosomes from the gtf files
  
  anno <- filter(anno_gtf, X1 %in% chroms, X3 == "transcript")
  anno$TXNAME <- gsub("^.+transcript_id \"(\\S+)\";.+;$", "\\1", anno$X9)
  anno <- unique.data.frame(select(anno, X1, X4, X5, X7, TXNAME))
  
  anno <- inner_join(tx2gene_local, anno, by = "TXNAME")
  # 7960 -> 6878 after this
  # QUESTION: will this cause loss of genes/transcripts during tximport?
  
  return(anno)
}

### apply the function
annoXY <- ChromAnno()
annoX <- ChromAnno(chroms = "chrX")
annoY <- ChromAnno(chroms = "chrY")



# predict the sample gender using X, Y gene expression --------------------

### function

OverviewXY <- function(txi, colData, brain_region,
                          design = ~ sex + condition,
                          pca_intergroup = "sex",
                          ggplot2_aes = aes(PC1, PC2, color = sex),
                          XYgenes, Xgenes, Ygenes){
  # modified fron the function DEanalysisAndSoOn()
  # use the expression of genes on chromosome X and T to predict the gender of the sample
  
  # ONLY for MSBB: remove samples from the colData, that are excluded from txi
  colData <- colData[rownames(colData) %in% colnames(txi$abundance), ]
  
  # remove samples from the txi, that are excluded from colData
  txi$abundance <- txi$abundance[, rownames(colData)]
  txi$counts <- txi$counts[, rownames(colData)]
  txi$length <- txi$length[, rownames(colData)]
  
  # ONLY for MSBB:
  # to avoid assay colnames() must be NULL or identical to colData rownames()
  colnames(txi$abundance) <- NULL
  colnames(txi$counts) <- NULL
  colnames(txi$length) <- NULL
  
  # subsetting txi by genes on chromosome X or Y
  txiX <- txi
  txiY <- txi
  
  # subset by gene on chromosome X
  txiX$abundance <- txi$abundance[row.names(txi$abundance) %in% Xgenes$GENEID, ]
  txiX$counts <- txi$counts[row.names(txi$counts) %in% Xgenes$GENEID, ]
  txiX$length <- txi$length[row.names(txi$length) %in% Xgenes$GENEID, ]
  
  # subset by gene on chromosome Y
  txiY$abundance <- txi$abundance[row.names(txi$abundance) %in% Ygenes$GENEID, ]
  txiY$counts <- txi$counts[row.names(txi$counts) %in% Ygenes$GENEID, ]
  txiY$length <- txi$length[row.names(txi$length) %in% Ygenes$GENEID, ]
  
  
  # construct dds and vsd object for genes on chromosome X
  ddsX <- DESeqDataSetFromTximport(txi = txiX, colData = colData,
                                   design = design)
  # collapse technical replicates
  ddsX <- collapseReplicates(object = ddsX, groupby = factor(colData[, "BioSample"]),
                             run = colData[, "Run"])
  vsdX <- varianceStabilizingTransformation(ddsX)
  
  # construct dds and vsd object for genes on chromosome Y
  ddsY <- DESeqDataSetFromTximport(txi = txiY, colData = colData,
                                  design = design)
  # collapse technical replicates
  ddsY <- collapseReplicates(object = ddsY, groupby = factor(colData[, "BioSample"]),
                             run = colData[, "Run"])
  vsdY <- varianceStabilizingTransformation(ddsY)
  
  
  # pca
  pca_dataY <- DESeq2::plotPCA(vsdY, intgroup = pca_intergroup,
                              returnData = TRUE)
  ggplot(data = pca_dataY, ggplot2_aes) +
    geom_point(size=3) + 
    geom_text(mapping = aes(label = name), size = 2) +
    ggsave(paste0("pca_", brain_region, "_chrY", ".png"))
  
  
  # calculate the ratio: total counts of genes on chromosome Y and that on chromosome X
  # for each sample (not each run)
  RatioYX <- apply(DESeq2::counts(ddsY), 2, sum) / apply(DESeq2::counts(ddsX), 2, sum)
  
  return(list(RatioYX = RatioYX, pca_dataY = pca_dataY))
}


DecideGender <- function(colData, RatioYX, pca_dataY, 
                         thresholdPC1 = mean(c(min(pca_dataY$PC1), max(pca_dataY$PC1)))){
  # after observing the PCA plot, 
  # manually choose a thresholds of the PC1 to define 2 groups;
  # determine which group represent which gender by comparing their mean RatioYX
  # added predicted sex to colData
  # NOTE: mean(n1, n2) returns n1, only mean(c(n1, n2)) return the average of n1, n2
  
  sex_infer <- matrix(nrow = dim(colData)[1], ncol = 1)
  colnames(sex_infer) <- "sex_infer"
  colData <- cbind(colData, sex_infer)
  
  avgRatioYX_a <- mean(RatioYX[pca_dataY$name[pca_dataY$PC1 < thresholdPC1]])
  avgRatioYX_b <- mean(RatioYX[pca_dataY$name[pca_dataY$PC1 > thresholdPC1]])
  print(paste("threshold for PC1 is ", thresholdPC1))
  print(paste("group on the left has an average Y - X ratio ", avgRatioYX_a))
  print(paste("group on the right has an average Y - X ratio ", avgRatioYX_b))
  
  if (avgRatioYX_a < avgRatioYX_b) {
    colData[colData[, "BioSample"] %in% pca_dataY$name[pca_dataY$PC1 < thresholdPC1], "sex_infer"] <- "F"
    colData[colData[, "BioSample"] %in% pca_dataY$name[pca_dataY$PC1 > thresholdPC1], "sex_infer"] <- "M"
    print("F on the left, M on the right")
  }
  else{
    colData[colData[, "BioSample"] %in% pca_dataY$name[pca_dataY$PC1 < thresholdPC1], "sex_infer"] <- "M"
    colData[colData[, "BioSample"] %in% pca_dataY$name[pca_dataY$PC1 > thresholdPC1], "sex_infer"] <- "F"
    print("M on the left, F on the right")
  }
  
  return(colData)
}


### apply the functions to different danasets

setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna")

# ----------- a tricky multi-line comment begin
if (FALSE) {


# GSE110731
txi_gse110731 <- readRDS("txi_GSE110731.rds")
colData_gse110731 <- readRDS("colData_GSE110731.rds")

# gender_pca_data_gse110731 <- 
OverviewXY_gse110731 <- OverviewXY(txi = txi_gse110731, colData = colData_gse110731,
                                brain_region = "GSE110731", 
                                design = ~ sex + age + condition,
                                XYgenes = annoXY, Xgenes = annoX, Ygenes = annoY)

colData_gse110731 <- DecideGender(colData = colData_gse110731,
                                  RatioYX = OverviewXY_gse110731$RatioYX,
                                  pca_dataY = OverviewXY_gse110731$pca_dataY)

saveRDS(colData_gse110731, "colData_GSE110731.rds")


# GSE53697
txi_gse53697 <- readRDS("txi_GSE53697.rds")
colData_gse53697 <- readRDS("colData_GSE53697.rds")

OverviewXY_gse53697 <- OverviewXY(txi = txi_gse53697, colData = colData_gse53697,
                                  brain_region = "GSE53697", design = ~ condition,
                                  ggplot2_aes = aes(PC1, PC2),
                                  XYgenes = annoXY, Xgenes = annoX, Ygenes = annoY)

colData_gse53697 <- DecideGender(colData = colData_gse53697,
                                 RatioYX = OverviewXY_gse53697$RatioYX,
                                 pca_dataY = OverviewXY_gse53697$pca_dataY)

saveRDS(colData_gse53697, "colData_GSE53697.rds")


# GSE125050 (must be devided into 4 cell types, low priority)
txi_gse125050 <- readRDS("txi_GSE125050.rds")
colData_gse125050 <- readRDS("colData_GSE125050.rds")

# OverviewXY_gse125050 <- OverviewXY(txi = txi_gse125050, colData = colData_gse125050,
#                                    brain_region = "GSE125050",
#                                    design = ~ sex + age + condition,
#                                    XYgenes = annoXY, Xgenes = annoX, Ygenes = annoY)
# 
# colData_gse125050 <- DecideGender(colData = colData_gse125050,
#                                  RatioYX = OverviewXY_gse125050$RatioYX,
#                                  pca_dataY = OverviewXY_gse125050$pca_dataY)

# saveRDS(colData_gse125050, "colData_GSE125050.rds")


# GSE104704
txi_gse104704 <- readRDS("txi_GSE104704.rds")
colData_gse104704 <- readRDS("colData_GSE104704.rds")

OverviewXY_gse104704 <- OverviewXY(txi = txi_gse104704, colData = colData_gse104704,
                                   brain_region = "GSE104704", design = ~ age + condition,
                                   ggplot2_aes = aes(PC1, PC2),
                                   XYgenes = annoXY, Xgenes = annoX, Ygenes = annoY)

colData_gse104704 <- DecideGender(colData = colData_gse104704,
                                  RatioYX = OverviewXY_gse104704$RatioYX,
                                  pca_dataY = OverviewXY_gse104704$pca_dataY)

saveRDS(colData_gse104704, "colData_GSE104704.rds")


# GSE95587
txi_gse95587 <- readRDS("txi_GSE95587.rds")
colData_gse95587 <- readRDS("colData_GSE95587.rds")

OverviewXY_gse95587 <- OverviewXY(txi = txi_gse95587, colData = colData_gse95587,
                                  brain_region = "GSE95587",
                                  design = ~ sex + age + condition,
                                  XYgenes = annoXY, Xgenes = annoX, Ygenes = annoY)

colData_gse95587 <- DecideGender(colData = colData_gse95587,
                                 RatioYX = OverviewXY_gse95587$RatioYX,
                                 pca_dataY = OverviewXY_gse95587$pca_dataY)

saveRDS(colData_gse95587, "colData_GSE95587.rds")


# GSE125583 (NOTE: SAMN10810551 is strange)
txi_gse125583 <- readRDS("txi_GSE125583.rds")
colData_gse125583 <- readRDS("colData_GSE125583.rds")

OverviewXY_gse125583 <- OverviewXY(txi = txi_gse125583, colData = colData_gse125583,
                                brain_region = "GSE125583",
                                design = ~ sex + age + condition,
                                XYgenes = annoXY, Xgenes = annoX, Ygenes = annoY)

colData_gse125583 <- DecideGender(colData = colData_gse125583,
                                 RatioYX = OverviewXY_gse125583$RatioYX,
                                 pca_dataY = OverviewXY_gse125583$pca_dataY)

saveRDS(colData_gse125583, "colData_GSE125583.rds")

# ROSMAP (FL)
txi_ROSMAP_FL <- readRDS("txi_ROSMAP.rds")
colData_ROSMAP_FL <- readRDS("colData_ROSMAP.rds")


OverviewXY_ROSMAP_FL <- OverviewXY(txi = txi_ROSMAP_FL, colData = colData_ROSMAP_FL,
                                   brain_region = "ROSMAP_FL",
                                   design = ~ sex + age + condition,
                                   XYgenes = annoXY, Xgenes = annoX, Ygenes = annoY)

colData_ROSMAP_FL <- DecideGender(colData = colData_ROSMAP_FL,
                                  RatioYX = OverviewXY_ROSMAP_FL$RatioYX,
                                  pca_dataY = OverviewXY_ROSMAP_FL$pca_dataY,
                                  thresholdPC1 = 5)

saveRDS(colData_ROSMAP_FL, "colData_ROSMAP_FL.rds")


# MSBB (FL)
txi_MSBB <- readRDS("txi_MSBB.rds")
colData_MSBB_FL <- readRDS("colData_MSBB_FL.rds")

OverviewXY_MSBB_FL <- OverviewXY(txi = txi_MSBB, colData = colData_MSBB_FL,
                                   brain_region = "MSBB_FL",
                                   design = ~ sex + age + condition,
                                   XYgenes = annoXY, Xgenes = annoX, Ygenes = annoY)

colData_MSBB_FL <- DecideGender(colData = colData_MSBB_FL,
                                  RatioYX = OverviewXY_MSBB_FL$RatioYX,
                                  pca_dataY = OverviewXY_MSBB_FL$pca_dataY)

colData_MSBB_FL <- colData_MSBB_FL[!is.na(colData_MSBB_FL[, 'sex_infer']), ]
# QUESTION: some samples passed patient inclusion criteria and are included by the colData,
# but not by the salmon folder. The reason is unknown.
saveRDS(colData_MSBB_FL, "colData_MSBB_FL.rds")


# MSBB (TL)
colData_MSBB_TL <- readRDS("colData_MSBB_TL.rds")


OverviewXY_MSBB_TL <- OverviewXY(txi = txi_MSBB, colData = colData_MSBB_TL,
                                   brain_region = "MSBB_TL",
                                   design = ~ sex + age + condition,
                                   XYgenes = annoXY, Xgenes = annoX, Ygenes = annoY)

colData_MSBB_TL <- DecideGender(colData = colData_MSBB_TL,
                                  RatioYX = OverviewXY_MSBB_TL$RatioYX,
                                  pca_dataY = OverviewXY_MSBB_TL$pca_dataY)

colData_MSBB_TL <- colData_MSBB_TL[!is.na(colData_MSBB_TL[, 'sex_infer']), ]
saveRDS(colData_MSBB_TL, "colData_MSBB_TL.rds")




# MayoRNAseq (CBE)

txi_MayoRNAseq_CBE <- readRDS("txi_MayoRNAseq_CB.rds")
colData_MayoRNAseq_CBE <- readRDS("colData_MayoRNAseq_CBE.rds")

OverviewXY_MayoRNAseq_CBE <- OverviewXY(txi = txi_MayoRNAseq_CBE, colData = colData_MayoRNAseq_CBE,
                                 brain_region = "MayoRNAseq_CBE",
                                 design = ~ sex + age + condition,
                                 XYgenes = annoXY, Xgenes = annoX, Ygenes = annoY)

colData_MayoRNAseq_CBE <- DecideGender(colData = colData_MayoRNAseq_CBE,
                                RatioYX = OverviewXY_MayoRNAseq_CBE$RatioYX,
                                pca_dataY = OverviewXY_MayoRNAseq_CBE$pca_dataY)

colData_MayoRNAseq_CBE <- colData_MayoRNAseq_CBE[!is.na(colData_MayoRNAseq_CBE[, 'sex_infer']), ]
# NOTE: 28 samples that passed the patient inclusion criteria and with broken fastq 
# are excluded temporarily

saveRDS(colData_MayoRNAseq_CBE, "colData_MayoRNAseq_CBE.rds")


}
# ----------- a tricky multi-line comment end



# MayoRNAseq (TCX)

txi_MayoRNAseq_TCX <- readRDS("txi_MayoRNAseq_TCX.rds")
colData_MayoRNAseq_TCX <- readRDS("colData_MayoRNAseq_TCX.rds")

OverviewXY_MayoRNAseq_TCX <- OverviewXY(txi = txi_MayoRNAseq_TCX, 
                                        colData = colData_MayoRNAseq_TCX,
                                        brain_region = "MayoRNAseq_TCX",
                                        design = ~ sex + age + condition,
                                        XYgenes = annoXY, Xgenes = annoX, Ygenes = annoY)

colData_MayoRNAseq_TCX <- DecideGender(colData = colData_MayoRNAseq_TCX,
                                       RatioYX = OverviewXY_MayoRNAseq_TCX$RatioYX,
                                       pca_dataY = OverviewXY_MayoRNAseq_TCX$pca_dataY)

colData_MayoRNAseq_TCX <- colData_MayoRNAseq_TCX[!is.na(colData_MayoRNAseq_TCX[, 'sex_infer']), ]
# no samples excluded

saveRDS(colData_MayoRNAseq_TCX, "colData_MayoRNAseq_TCX.rds")
