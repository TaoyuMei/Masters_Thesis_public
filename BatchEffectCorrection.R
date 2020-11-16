# Batch Effect Correction using different algorithms
# this script only include funciton
# to be sourced by DESeq2_EachBrainRegions.R

library(sva)
library(RUVSeq)


# correction of batch effect ----------------------------------------------

### covariates in the design formula
NumeriseCovariates <- function(dds, design){
  # colData matrix coerce numerical variables (i.e. covariates rather than batches)
  # into characters and then factors, causing errors in the model matrix,
  # which need to be reverse.
  
  col_vars <- colnames(colData(dds))
  if("age" %in% col_vars){  # convert from factor to real number, like SV1
    dds$age <- as.numeric(as.character(dds$age))
  }
  if ("pmi" %in% col_vars) {
    dds$pmi <- as.numeric(as.character(dds$pmi))
  }
  
  # more covariates to be added
  
  design(dds) <- design
  
  dds <- dds[, !is.na(dds$pmi)] 
  dds <- dds[, !is.na(dds$age)]
  
  return(dds)
}


### sva's svaseq
SVAtoDDS <- function(dds, n_sv){
  # apply svaseq to the dds object and return an adjusted one
  
  dds <- estimateSizeFactors(dds)  # prepare for normalisation

  dat  <- counts(dds, normalized = TRUE)
  
  idx  <- rowMeans(dat) > 1  #  errors will be caused without this
  dat  <- dat[idx, ]
  
  mod  <- model.matrix(~ condition, colData(dds))
  mod0 <- model.matrix(~ 1, colData(dds))
  svseq <- svaseq(dat, mod, mod0, n.sv = n_sv)
  
  
  if (n_sv == 1) {
    
    dds$SV1 <- svseq$sv[, 1]
    
    design(dds) <- ~ SV1 + condition
  }
  else if (n_sv == 2) {
    dds$SV1 <- svseq$sv[, 1]
    dds$SV2 <- svseq$sv[, 2]
    
    design(dds) <- ~ SV1 + SV2 + condition
  }
  else if (n_sv == 3) {
    dds$SV1 <- svseq$sv[, 1]
    dds$SV2 <- svseq$sv[, 2]
    dds$SV3 <- svseq$sv[, 3]
    
    design(dds) <- ~ SV1 + SV2 + SV3 + condition
  }
  else if (n_sv == 4) {
    dds$SV1 <- svseq$sv[, 1]
    dds$SV2 <- svseq$sv[, 2]
    dds$SV3 <- svseq$sv[, 3]
    dds$SV4 <- svseq$sv[, 4]
    
    design(dds) <- ~ SV1 + SV2 + SV3 + SV4 + condition
  }
  
  return(dds)
}


BatchesSVAtoDDS <- function(dds, n_sv){
  # apply svaseq to the dds object and return an adjusted one
  
  dds <- estimateSizeFactors(dds)  
  dat  <- counts(dds, normalized = TRUE)
  
  idx  <- rowMeans(dat) > 1  
  dat  <- dat[idx, ]
  
  mod  <- model.matrix(design(dds), colData(dds))  # e.g. ~ age + pmi + condition
  mod0 <- model.matrix(~ 1, colData(dds))
  
  svseq <- svaseq(dat, mod, mod0, n.sv = n_sv)  
  
  if (n_sv == 1) {
    dds$SV1 <- svseq$sv[, 1]
    
    for_ele <- strsplit(toString(design(dds)), split = "[,|+]")
    nele <- length(for_ele[[1]])
    design(dds) <- as.formula(paste(c(paste(for_ele[[1]][1:2], collapse = ""),
                                      for_ele[[1]][3:(nele - 1)], 
                                      "SV1", 
                                      for_ele[[1]][nele]), collapse = "+"))
  }
  else if (n_sv == 2) {
    dds$SV1 <- svseq$sv[, 1]
    dds$SV2 <- svseq$sv[, 2]
    
    for_ele <- strsplit(toString(design(dds)), split = "[,|+]")
    nele <- length(for_ele[[1]])
    design(dds) <- as.formula(paste(c(paste(for_ele[[1]][1:2], collapse = ""),
                                      for_ele[[1]][3:(nele - 1)], 
                                      "SV1", "SV2",
                                      for_ele[[1]][nele]), collapse = "+"))
  }
  else if (n_sv == 3) {
    dds$SV1 <- svseq$sv[, 1]
    dds$SV2 <- svseq$sv[, 2]
    dds$SV3 <- svseq$sv[, 3]
    
    for_ele <- strsplit(toString(design(dds)), split = "[,|+]")
    nele <- length(for_ele[[1]])
    design(dds) <- as.formula(paste(c(paste(for_ele[[1]][1:2], collapse = ""),
                                      for_ele[[1]][3:(nele - 1)], 
                                      "SV1", "SV2", "SV3",
                                      for_ele[[1]][nele]), collapse = "+"))
  }
  else if (n_sv == 4) {
    dds$SV1 <- svseq$sv[, 1]
    dds$SV2 <- svseq$sv[, 2]
    dds$SV3 <- svseq$sv[, 3]
    dds$SV4 <- svseq$sv[, 4]
    
    for_ele <- strsplit(toString(design(dds)), split = "[,|+]")
    nele <- length(for_ele[[1]])
    design(dds) <- as.formula(paste(c(paste(for_ele[[1]][1:2], collapse = ""),
                                      for_ele[[1]][3:(nele - 1)], 
                                      "SV1", "SV2", "SV3", "SV4",
                                      for_ele[[1]][nele]), collapse = "+"))
  }
  
  return(dds)
}




### RUVseq
NegativeControlGenesFromRes <- function(ddsM, resM, ddsF, resF){
  # get negative control genes from DESeq results
  
  setM <- newSeqExpressionSet(counts(ddsM))
  idxM  <- rowSums(counts(setM) > 5) >= 2
  setM  <- setM[idxM, ]
  setM <- betweenLaneNormalization(setM, which = "upper")
  not.sigM <- rownames(resM)[which(resM$pvalue > 0.1)]
  empiricalM <- rownames(setM)[rownames(setM) %in% not.sigM]
  
  setF <- newSeqExpressionSet(counts(ddsF))
  idxF  <- rowSums(counts(setF) > 5) >= 2
  setF  <- setF[idxF, ]
  setF <- betweenLaneNormalization(setF, which = "upper")
  not.sigF <- rownames(resF)[which(resF$pvalue > 0.1)]
  empiricalF <- rownames(setF)[rownames(setF) %in% not.sigF]
  
  return(list(empiricalM = empiricalM, empiricalF = empiricalF))
}

GetWsByRUVg <- function(empiricalM, empiricalF,
                        ddsM, resM, ddsF, resF){
  
  setM <- newSeqExpressionSet(counts(ddsM))
  idxM  <- rowSums(counts(setM) > 5) >= 2
  setM  <- setM[idxM, ]
  setM <- betweenLaneNormalization(setM, which = "upper")
  
  setF <- newSeqExpressionSet(counts(ddsF))
  idxF  <- rowSums(counts(setF) > 5) >= 2
  setF  <- setF[idxF, ]
  setF <- betweenLaneNormalization(setF, which = "upper")
  
  setM <- RUVg(setM, empiricalM, k=2)
  setF <- RUVg(setF, empiricalF, k=2)
  
  return(list(M_W1 = setM$W_1, M_W2 = setM$W_2, 
              F_W1 = setF$W_1, F_W2 = setF$W_2))
}
