#!/usr/local/bin/Rscript-3.6.2
# preparing for DE analysis, using tximport
library(tximport)
library(tidyverse)



# prepare to apply tximport to the complete datasets  -------------------------------

## prepare the transcripts - gene table

setwd("/binf-isilon/alab/students/vrw936/scratch/reference_ome_and_anno/")
gencode.v33.metadata.HGNC <- read_tsv("gencode.v33.metadata.HGNC",
                                      col_names = FALSE)
gencode.v33.metadata.EntrezGene <- read_tsv("gencode.v33.metadata.EntrezGene",
                                            col_names = FALSE)

tx2gene <- select(gencode.v33.metadata.HGNC, 
                  TXNAME = X1,
                  GENEID = X2)

setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna")


## functionalise

txiForEachBrainRegion <- function(tx2gene, 
                                  wd = "/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna", 
                                  brainRegion, datasets_dir, datasetType = "GSE"){
  # a function applying tximport to all salmon quant results of a given brain region
  # tx2gene: a table containing the IDs of all transcripts 
  # that are included in the reference transcriptome, 
  # and the ID/Symbol of the gene that each transcript belongs to.
  # wd: working directory containing the datasets' folders.
  # brainRegion: a character string describing the brain region, without space
  # datasets_dir: relevent directories of all folders containing fastq files
  
  setwd(wd)
  if(datasetType == "GSE"){
    srrFolders <- list.files(file.path(".", datasets_dir, "trimmed", "salmon"),
                             full.names = TRUE)
  }
  else if(datasetType == "AMP_AD"){
    srrFolders <- list.files(file.path(".", datasets_dir, "trimmed", "shuffled", "salmon"),
                             full.names = TRUE)
  }
  
  # AMP-AD
  
  
  quantFiles <- file.path(srrFolders, "quant.sf")
  txi <- tximport(quantFiles, type = "salmon", tx2gene = tx2gene)  
  
  # GSE
  # srrs <- list.files(file.path(".", datasets_dir, 
  #                              "trimmed", "salmon"))
  # AMP-AD
  srrs <- list.files(file.path(".", datasets_dir, "trimmed", "shuffled", "salmon"))
  
  colnames(txi$counts) <- srrs
  colnames(txi$abundance) <- srrs
  colnames(txi$length) <- srrs
  
  saveRDS(txi, paste0("txi_", brainRegion, ".rds"))
  
  # return(txi)
}



# apply tximport within each dataset --------------------------------------

txiForEachBrainRegion(tx2gene = tx2gene, brainRegion = "GSE104704",
                      datasets_dir = "GSE104704")

txiForEachBrainRegion(tx2gene = tx2gene, brainRegion = "GSE110731",
                      datasets_dir = "GSE110731")

txiForEachBrainRegion(tx2gene = tx2gene, brainRegion = "GSE125050",
                      datasets_dir = "GSE125050")

txiForEachBrainRegion(tx2gene = tx2gene, brainRegion = "GSE125583",
                      datasets_dir = "GSE125583")

txiForEachBrainRegion(tx2gene = tx2gene, brainRegion = "GSE95587",
                      datasets_dir = "GSE95587")

txiForEachBrainRegion(tx2gene = tx2gene, brainRegion = "GSE53697",
                      datasets_dir = "GSE53697")

txiForEachBrainRegion(tx2gene = tx2gene, brainRegion = "ROSMAP",
                      datasets_dir = "ROSMAP/fastq")

txiForEachBrainRegion(tx2gene = tx2gene, brainRegion = "MayoRNAseq_TCX",
                      datasets_dir = "MayoRNAseq/fastq/Mayo_TCX_sample_FASTQs/")

txiForEachBrainRegion(tx2gene = tx2gene, brainRegion = "MayoRNAseq_CB",
                      datasets_dir = "MayoRNAseq/fastq/Mayo_CBE_sample_FASTQs/")

txiForEachBrainRegion(tx2gene = tx2gene, brainRegion = "MSBB",
                      datasets_dir = "MSBB/fastq")
# NOTE: one txi for MSBB, two colDatas for TL and FL respectively
