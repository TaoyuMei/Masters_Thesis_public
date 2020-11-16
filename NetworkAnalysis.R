# network analysis
library(limma)
library(tidyverse)
library(RCy3)


# functions ---------------------------------------------------------------



STRING_visual_NewMembers <- function(fig_dir = "STRING",
                                     DE_result, brain_region, 
                                     lfc_threshold = 0.1, alpha = 0.05,
                                     fig_name, key_genes){
  
  
  dir.create(fig_dir)
  
  DE_result <- filter(as_tibble(DE_result, rownames = "gene"), 
                      padj < alpha, abs(log2FoldChange) > lfc_threshold)
  
  key_genes <- key_genes[key_genes %in% DE_result$gene]
  
  # split gene list into chunks due to the limitation on url length of the STRING server
  max <- 600 - length(key_genes)
  candid_new_mem <- DE_result$gene[!(DE_result$gene %in% key_genes)]
  x <- seq_along(candid_new_mem)
  candid_new_mem <- split(candid_new_mem, ceiling(x/max))
  
  key_genes_string <- paste0(key_genes, collapse = "%0d")
  n = 1
  interactors <- c()
  
  for(chunk in candid_new_mem){
    chunk <- paste0(chunk, collapse = "%0d")
    your_identifiers <- paste0(chunk, "%0d", key_genes_string)
    
    STRING_ppi_url <- paste0("https://string-db.org/api/",
                             "tsv",
                             "/network?identifiers=",
                             your_identifiers,
                             "&species=9606")
    
    download.file(STRING_ppi_url, 
                  destfile = file.path(fig_dir, paste0(brain_region, "_", fig_name,
                                                       "_", n, "_ppi.tsv")))
    
    
    tsv_tmp <- read_tsv(file.path(fig_dir, paste0(brain_region, "_", 
                                                  fig_name, "_", n, "_ppi.tsv")))
    
    A_B_tmp <- c()
    for (i in 1:dim(tsv_tmp)[1]) {
      if((tsv_tmp$preferredName_A[i] %in% key_genes) &
         !(tsv_tmp$preferredName_B[i] %in% key_genes)){
        
        A_B_tmp <- c(A_B_tmp, tsv_tmp$preferredName_B[i])
      }
      if((tsv_tmp$preferredName_B[i] %in% key_genes) &
         !(tsv_tmp$preferredName_A[i] %in% key_genes)){
        
        A_B_tmp <- c(A_B_tmp, tsv_tmp$preferredName_A[i])
      }
    }
    
    interactors <- c(interactors, A_B_tmp)
    
    n <- n + 1
    
  }
  
  # table(table(interactors))
  # new members should interact with at least one third of the old members
  
  threshold_number <- round(length(key_genes) / 3)  
  genes_selected <- names(table(interactors)[table(interactors) >= threshold_number])
  
  # genes_selected <- names(sort(table(interactors), decreasing = TRUE)[1:20])
  
  
  system(paste0("rm ", file.path(fig_dir, 
                                 paste0(brain_region, "_", fig_name, 
                                        "_", "*", "_ppi.tsv"))))
  
  rm(list = ls(pattern = paste0(brain_region, "_", fig_name,
                                "_", ".+", "_ppi")))
  
  genes_selected_str <- paste0(genes_selected, collapse = "%0d")
  your_identifiers <- paste0(genes_selected_str, "%0d", key_genes_string)
  
  STRING_ppi_url <- paste0("https://string-db.org/api/",
                           "tsv",
                           "/network?identifiers=",
                           your_identifiers,
                           "&species=9606")
  
  download.file(STRING_ppi_url, 
                destfile = file.path(fig_dir, paste0(brain_region, "_", 
                                                     fig_name, "_NewMember_ppi.tsv")))
  
  # e.g. combined_FL_M_mitophagy_NewMember_ppi
  
  ppi_tmp <- read_tsv(file.path(fig_dir,
                                paste0(brain_region, "_", fig_name, "_NewMember_ppi.tsv")))
  
  ppi_tmp <- dplyr::filter(ppi_tmp, 
                           preferredName_A %in% DE_result$gene, 
                           preferredName_B %in% DE_result$gene)
  
  # add the degree as an attribute of each node (gene) in 2 new column
  # add the log2FC, and whether it is a new member  
  ppi_tmp$A_lfc <- 0
  ppi_tmp$B_lfc <- 0
  ppi_tmp$A_degree <- 0
  ppi_tmp$B_degree <- 0
  ppi_tmp$whether_Anew <- "old"
  ppi_tmp$whether_Bnew <- "old"
  
  node_degree <- table(c(ppi_tmp$preferredName_A, ppi_tmp$preferredName_B))
  
  for (j in 1:dim(ppi_tmp)[1]) {
    ppi_tmp$A_degree[j] <- as.numeric(node_degree[ppi_tmp$preferredName_A[j]])
    ppi_tmp$A_lfc[j] <- DE_result[DE_result$gene == ppi_tmp$preferredName_A[j], 
                                  "log2FoldChange"][[1]]
    if(ppi_tmp$preferredName_A[j] %in% genes_selected){
      ppi_tmp$whether_Anew[j] <- "new"
    }
    
    ppi_tmp$B_degree[j] <- as.numeric(node_degree[ppi_tmp$preferredName_B[j]])
    ppi_tmp$B_lfc[j] <- DE_result[DE_result$gene == ppi_tmp$preferredName_B[j], 
                                  "log2FoldChange"][[1]]
    if(ppi_tmp$preferredName_B[j] %in% genes_selected){
      ppi_tmp$whether_Bnew[j] <- "new"
    }
  }
  
  ppi_tmp <- dplyr::select(ppi_tmp, preferredName_A, preferredName_B, 
                           A_lfc, B_lfc, A_degree, B_degree, whether_Anew, whether_Bnew)
  
  write_tsv(x = ppi_tmp, path = file.path(fig_dir,
                                          paste0(brain_region, "_", fig_name, "_NewMember_ppi.tsv")))
  
  
  # also download the 'oficial' figure
  STRING_image_url <- paste0("https://string-db.org/api/",
                             "image",
                             "/network?identifiers=",
                             your_identifiers,
                             "&species=9606",
                             "&hide_disconnected_nodes=1")
  
  download.file(STRING_image_url, 
                destfile = file.path(fig_dir, paste0(brain_region, 
                                                     "_", fig_name, "_NewMember_ppi.png")))
  
}



# test the function
# 
STRING_visual_NewMembers(DE_result = DEres_combined_FL_M_braak_batches_sva_vst,
                         brain_region = "combined_FL_M",
                         fig_name = "mitophagy",
                         key_genes = c(combined_FL_M_down_mitophagy_DEGs,
                                       combined_FL_M_up_mitophagy_DEGs))



# apply the functions to each datasets/brain regions -----------------------


##### STRING

setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna")

# load DEres objects and key gene lists first
dds_combined_FL_F_batches_sva_vst <- readRDS("DE_combined_batches_sva_figures/dds_combined_FL_F_batches_sva_vst.rds")
dds_combined_FL_M_batches_sva_vst <- readRDS("DE_combined_batches_sva_figures/dds_combined_FL_M_batches_sva_vst.rds")
DEres_combined_FL_F_braak_batches_sva_vst <- readRDS("DE_combined_batches_sva_figures/DEres_combined_FL_F_braak_batches_sva_vst.rds")
DEres_combined_FL_M_braak_batches_sva_vst <- readRDS("DE_combined_batches_sva_figures/DEres_combined_FL_M_braak_batches_sva_vst.rds")


geneset_dir = "GeneSet_figures"

for (deg in list.files(path = geneset_dir, pattern = "_DEGs")) {
  assign(gsub("(^.+)\\.rds", "\\1", deg), 
         readRDS(file.path(geneset_dir, deg)))
}

# e.g. combined_FL_F_down_lysosome_DEGs

for (deg in ls(pattern = "combined.+_DEGs")){  # M and F are separated
  sex_info <- gsub("combined_.+_([F|M])_.+_.+_DEGs", "\\1", deg)
  dataset_info <- gsub("(combined_.+_[F|M])_.+_.+_DEGs", "\\1", deg)
  up_or_down <- gsub("combined_.+_[F|M]_(.+)_.+_DEGs", "\\1", deg)
  fig_name <- gsub("combined_.+_[F|M]_.+_(.+)_DEGs", "\\1", deg)
  
  if (up_or_down == "up") {
    
    key_genes <- c(get(deg), get(paste0(dataset_info, "_down_", fig_name, "_DEGs")))
    # combine up- and down- regulated genes
    
    STRING_visual_NewMembers(DE_result = get(paste0("DEres_", dataset_info, 
                                                    "_braak_batches_sva_vst")),
                             brain_region = dataset_info,
                             fig_name = fig_name,
                             key_genes = key_genes)
  }
  
}



# extract and count the highly-connected genes ----------------------------


for(ppi_tsv in list.files(path = "STRING", pattern = "tsv")){
  sex_info <- gsub("combined_.+_([F|M])_.+_.+_ppi.tsv", "\\1", ppi_tsv)
  
  if (sex_info == "M") {
    print("")
    print("")
    print(ppi_tsv)
    
    # M
    tsv_tmp <- read_tsv(file.path("STRING", ppi_tsv))
    
    A <- dplyr::filter(tsv_tmp, A_degree >= 10)$preferredName_A
    B <- dplyr::filter(tsv_tmp, B_degree >= 10)$preferredName_B
    all_nodes <- c(A, B)
    
    highly_cont_gene <- unique.default(all_nodes)
    
    print(sex_info)
    print(length(highly_cont_gene))
    print(paste0(highly_cont_gene, collapse = ", "))
    
    # F
    ppi_tsv_F <- sub("_M_", "_F_", ppi_tsv)
    tsv_tmp_F <- read_tsv(file.path("STRING", ppi_tsv_F))
    
    A <- dplyr::filter(tsv_tmp_F, A_degree >= 10)$preferredName_A
    B <- dplyr::filter(tsv_tmp_F, B_degree >= 10)$preferredName_B
    all_nodes <- c(A, B)
    
    highly_cont_gene_F <- unique.default(all_nodes)
    
    print("F")
    print(length(highly_cont_gene_F))
    print(paste0(highly_cont_gene_F, collapse = ", "))
    
    
    # M & F
    overlap <- highly_cont_gene[highly_cont_gene %in% highly_cont_gene_F]
    print("overlap")
    print(length(overlap))
    print(paste0(overlap, collapse = ", "))
  }
}





