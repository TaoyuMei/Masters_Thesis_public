# put together all functions related to the use of DESeq2

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(tidyverse)
library(VennDiagram)
source('/binf-isilon/alab/students/vrw936/Master_Thesis/MyCode/BatchEffectCorrection.R')


# 2 functions to construct DDS objects for 2 sexes, do DE analysis and visualise the result -----------


ConstructDDSandVisualise <- function(txi, colData, brain_region, fig_dir = "./DE_sva_figures/",
                                     BEC = FALSE, correct_condition = TRUE, transform = "vst",
                                     tree_or_not = TRUE, Ws = NULL, n_sv = NULL,
                                     design = ~ condition, pca_intergroup = c("condition"),
                                     ggplot2_aes = aes(PC1, PC2, color = condition)){
  # to detect & correct batch effect
  
  
  
  dir.create(fig_dir)  # no problem even if the folder exists
  
  # correct the condition according to braak score 
  # only applicable to GSE125583, GSE125050, GSE95587, MayoRNAseq
  # should also affect the file name saved
  figureSuffix <- ""
  if (correct_condition == TRUE) {
    colData[colData[, "braak"] %in% c("IV", "V", "VI", "4.5", "5.0", "5.5", "6.0"), "condition"] <- "AD"
    colData[colData[, "braak"] %in% c("I", "II", "0.0", "1.0", "2.0", "2.5", NA), "condition"] <- "ctrl"
    colData[colData[, "braak"] %in% c("III", "3.0"), "condition"] <- "out"
    colData <- colData[colData[, "condition"] != "out", ]
    figureSuffix <- paste0(figureSuffix, "_braak")  
  }
  
  # remove samples from the colData, that are excluded from txi
  colData <- colData[rownames(colData) %in% colnames(txi$abundance), ]
  
  # remove samples from the txi, that are excluded from colData
  colData <- colData[order(colData[, "condition"], decreasing = TRUE), ]
  txi$abundance <- txi$abundance[, rownames(colData)]
  txi$counts <- txi$counts[, rownames(colData)]
  txi$length <- txi$length[, rownames(colData)]
  
  # to avoid "assay colnames() must be NULL or identical to colData rownames()"
  all(colnames(txi$abundance) == rownames(colData))
  all(colnames(txi$length) == rownames(colData))
  all(colnames(txi$counts) == rownames(colData))
  colnames(txi$abundance) <- NULL
  colnames(txi$counts) <- NULL
  colnames(txi$length) <- NULL
  
  
  # construction of the dds object
  # using counts and average transcript lengths from tximport
  dds <- DESeqDataSetFromTximport(txi = txi, colData = colData,
                                  design = ~ condition)

  # collapse technical replicates 
  dds <- collapseReplicates(object = dds, groupby = factor(colData[, "BioSample"]),
                            run = colData[, "Run"])
  
  # set reference level
  dds$condition <- relevel(dds$condition, ref = "ctrl")
  
  # divide male and female (remember that samples correspond to columns)
  ddsM <- dds[, dds$sex_infer == "M"]
  ddsF <- dds[, dds$sex_infer == "F"]
  
  ## Batch Effect Correction (call another function)
  if(BEC != FALSE){
    
    if (BEC == "sva") {
      ddsM <- SVAtoDDS(ddsM, n_sv)
      ddsF <- SVAtoDDS(ddsF, n_sv)
    }
    else if (BEC == "batches"){
      dds <- NumeriseCovariates(dds, design)
      ddsM <- NumeriseCovariates(ddsM, design)
      ddsF <- NumeriseCovariates(ddsF, design)
    }
    else if (BEC == "batches_sva") {
      dds <- NumeriseCovariates(dds, design)
      ddsM <- NumeriseCovariates(ddsM, design)
      ddsF <- NumeriseCovariates(ddsF, design)
      
      ddsM <- BatchesSVAtoDDS(ddsM, n_sv)
      ddsF <- BatchesSVAtoDDS(ddsF, n_sv)
    }
    else if (BEC == "RUVseq") {
      ddsM$W1 <- Ws$M_W1
      ddsF$W1 <- Ws$F_W1
      
      ddsM$W2 <- Ws$M_W2
      ddsF$W2 <- Ws$F_W2
      
      design(ddsM) <- ~ W1 + W2 + condition
      design(ddsF) <- ~ W1 + W2 + condition
    }
    else if (BEC == "comBat-seq") {
      
    }
    
    figureSuffix <- paste0(figureSuffix, "_", BEC)
  }
  else{  # BEC == FALSE
    transdM <- ddsM
    transdF <- ddsF
  }
  
  
  ## overview of the gene expression before DE analysis
  # using 'avgTxLength' from assays(dds), correcting for library size
  # transform the counts matrix before visualisation: choose between vst, rlt, nt
  if (transform == "vst") {
    
    transdM <-DESeq2::vst(ddsM)  # default: blind = TRUE
    modM <- model.matrix(design(ddsM), colData(ddsM))
    matM <- assay(transdM)
    
    transdF <-DESeq2::vst(ddsF)  # default: blind = TRUE
    modF  <- model.matrix(design(ddsF), colData(ddsF))
    matF <- assay(transdF)
    
    if (BEC %in% c("batches", "sva", "RUVseq", "batches_sva")) {
      matM <- limma::removeBatchEffect(x = matM,
                                       #design = model.matrix(~ condition, colData(transdM)),
                                       covariates = modM[, c(-1, -dim(modM)[2])])
      
      matF <- limma::removeBatchEffect(x = matF,
                                       #design = model.matrix(~ condition, colData(transdF)),
                                       covariates = modF[, c(-1, -dim(modF)[2])])
    
      print(paste(c("covariates are", paste0(colnames(modM)[c(-1, -dim(modM)[2])]), 
                                          collapse = ", ")))
    }
    
    assay(transdM) <- matM
    assay(transdF) <- matF
    
    figureSuffix <- paste0(figureSuffix, "_vst")
  }
  else if (transform == "norm") {  # normTransform() is not used in this way 
    transdM <- normTransform(ddsM)  
    transdF <- normTransform(ddsF)
    figureSuffix <- paste0(figureSuffix, "_norm")
  }
  else if (transform == "rlog") {
    transdM <- rlog(ddsM)
    transdF <- rlog(ddsF)
    figureSuffix <- paste0(figureSuffix, "_rlog")
  }
  
  # PCA  
  pca_dataM <- DESeq2::plotPCA(transdM, intgroup = pca_intergroup,
                               returnData = TRUE)
  ggplot(data = pca_dataM, ggplot2_aes) +
    geom_point(size=3) + 
    geom_text(mapping = aes(label = name), size = 2) +
    # labs(title = paste0("pca_", brain_region, "_M", figureSuffix)) +
    ggsave(paste0(fig_dir, "pca_", brain_region, "_M", figureSuffix, ".png"))
  
  pca_dataF <- DESeq2::plotPCA(transdF, intgroup = pca_intergroup,
                               returnData = TRUE)
  ggplot(data = pca_dataF, ggplot2_aes) +
    geom_point(size=3) + 
    geom_text(mapping = aes(label = name), size = 2) +
    ggsave(paste0(fig_dir, "pca_", brain_region, "_F", figureSuffix,".png"))
  
  
  # hierarchical clustering tree
  # it seems that assay(vsd) can get the expression matrix (not DESeq2::assay(), just assay())
  # assay() is also applicable to dds, the same as counts()
  if (tree_or_not == TRUE) {
    pdf(file = paste0(fig_dir, "tree_", brain_region, "_M", figureSuffix,".pdf"))
    expr.matr <- assay(transdM)
    expr.matr.t <- t(expr.matr)
    eu.distM <- dist(expr.matr.t, method = "euclidean")
    tree <- hclust(eu.distM, method = "average")
    plot(tree, main = "", xlab = "Sample", sub = "")
    dev.off()
    
    pdf(file = paste0(fig_dir, "tree_", brain_region, "_F", figureSuffix,".pdf"))
    expr.matr <- assay(transdF)
    expr.matr.t <- t(expr.matr)
    eu.distF <- dist(expr.matr.t, method = "euclidean")
    tree <- hclust(eu.distF, method = "average")
    plot(tree, main = "", xlab = "Sample", sub = "")
    dev.off()
  }
  
  
  # heatmap (select the top 100 genes that contribute most to the variance)
  pdf(file = paste0(fig_dir, "heatmap_", brain_region, "_M", figureSuffix,".pdf"))
  select <- order(rowMeans(assay(transdM)), decreasing = TRUE)[1:50]
  df <- as.data.frame(colData(ddsM)[, c("BioSample", pca_intergroup)])
  df <-  dplyr::select(df, condition)
  pheatmap(assay(transdM)[select,], cluster_rows = FALSE, show_rownames = FALSE,
           cluster_cols = TRUE, annotation_col = df,
           cutree_cols = 2)
  dev.off()
  
  pdf(file = paste0(fig_dir, "heatmap_", brain_region, "_F", figureSuffix,".pdf"))
  select <- order(rowMeans(assay(transdF)), decreasing = TRUE)[1:50]
  df <- as.data.frame(colData(ddsF)[, c("BioSample", pca_intergroup)])
  df <-  dplyr::select(df, condition)
  pheatmap(assay(transdF)[select,], cluster_rows = FALSE, show_rownames = FALSE,
           cluster_cols = TRUE, annotation_col = df, 
           cutree_cols = 2)
  dev.off()
  
  saveRDS(ddsM, paste0(fig_dir, "dds_", brain_region, "_M", figureSuffix, ".rds"))
  saveRDS(ddsF, paste0(fig_dir, "dds_", brain_region, "_F", figureSuffix, ".rds"))
  return(list(ddsM = ddsM, ddsF = ddsF))
}


DEanalysis <- function(ddsM, ddsF, brain_region, alpha = 0.05,
                       fig_dir = "./DE_sva_figures/", figureSuffix = ""){
  # differential analysis using DESeq2; 
  # use 2 batch effect-corrected DDSs as input, representing M and F
  # visualise the results
  
  system(paste("echo", brain_region, ": ", dim(colData(ddsM))[1], " male and ", 
               dim(colData(ddsF))[1], " female samples begin"))
  system("date")
  dir.create(fig_dir)  # no problem even if the folder exists
  
  
  ddsM <- DESeq(ddsM)
  ddsF <- DESeq(ddsF)
  
  resM <- results(ddsM, contrast = c('condition', 'AD', 'ctrl'))
  resF <- results(ddsF, contrast = c('condition', 'AD', 'ctrl'))

  # save the results
  saveRDS(resM, paste0(fig_dir, "DEres_", brain_region, "_M", figureSuffix, ".rds"))
  saveRDS(resF, paste0(fig_dir, "DEres_", brain_region, "_F", figureSuffix, ".rds"))
  
  
  ## visualise the results as volcano plots and MA plots
  
  # MA plot
  pdf(file = paste0(fig_dir, "MA_", brain_region, "_M", figureSuffix, ".pdf"))
  DESeq2::plotMA(resM, alpha = alpha, main = "")
  dev.off()
  
  pdf(file = paste0(fig_dir, "MA_", brain_region, "_F", figureSuffix, ".pdf"))
  DESeq2::plotMA(resF, alpha = alpha, main = "")
  dev.off()
  
  # volcano plot
  dev.new()
  pdf(file = paste0(fig_dir, "volano_", brain_region, "_M", figureSuffix, ".pdf"))
  EnhancedVolcano(resM, lab = rownames(resM), x = 'log2FoldChange',
                  y = 'pvalue', xlim = c(-10, 10), pCutoff = alpha, FCcutoff = 0.15)
  dev.off()
  
  dev.new()
  pdf(file = paste0(fig_dir, "volano_", brain_region, "_F", figureSuffix, ".pdf"))
  EnhancedVolcano(resF, lab = rownames(resF), x = 'log2FoldChange',
                  y = 'pvalue', xlim = c(-10, 10), pCutoff = alpha, FCcutoff = 0.15)
  dev.off()
  
  
  system(paste("echo", brain_region, ": ", dim(colData(ddsM))[1], " male and ", 
               dim(colData(ddsF))[1], " female samples finish"))
  system("date")
  
  # show the number of DEGs
  # print(paste0("Male: ", table(as_tibble(resM$padj)$value < 0.05)["TRUE"], " DEGs"))
  # print(paste0("Female: ", table(as_tibble(resF$padj)$value < 0.05)["TRUE"], " DEGs"))
  

}

