# make combined txi and colData for each brain regions, construct dds
# and finish following analyses
source('/binf-isilon/alab/students/vrw936/Master_Thesis/MyCode/DESeq2_functions.R')


# functions ---------------------------------------------------------------


Combine_ColDatas_Txis <- function(colDataList, txiList, brain_region, 
                                  col_vars = c("BioSample", "Run", "sex_infer", 
                                               "batch", "age", "pmi", "GSE", "condition")){
  # input: a list of dds objects, and a list of corresponding txi objects 
  # for the same brain region;
  # output: a combined dds object for that brain region, 
  # with batch effect varables properly described in colData and design formula;
  colData_combined <- colDataList[[1]][, col_vars]
  txi_combined <- txiList[[1]]
  
  for(i in 2:length(colDataList)){
    colData_combined <- rbind(colData_combined, 
                              colDataList[[i]][, col_vars])
    txi_combined$counts <- cbind(txi_combined$counts, txiList[[i]]$counts)
    txi_combined$length <- cbind(txi_combined$length, txiList[[i]]$length)
    txi_combined$abundance <- cbind(txi_combined$abundance, txiList[[i]]$abundance)
  }
  
  saveRDS(txi_combined, paste0("txi_combined_", brain_region, ".rds"))
  saveRDS(colData_combined, paste0("colData_combined_", brain_region, ".rds"))
  
  return(list(colData_combined = colData_combined, 
              txi_combined = txi_combined))
}




# apply functions to data -------------------------------------------------

setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna")


##### frontal lobe

colData_GSE110731 <- readRDS("colData_GSE110731.rds")
colData_MSBB_FL <- readRDS("colData_MSBB_FL.rds")
colData_ROSMAP_FL <- readRDS("colData_ROSMAP_FL.rds")

txi_GSE110731 <- readRDS("txi_GSE110731.rds")
txi_ROSMAP <- readRDS("txi_ROSMAP.rds")
txi_MSBB <- readRDS("txi_MSBB.rds")

# redefine batches according to PCA plots without adjustment
for (i in 1:dim(colData_MSBB_FL)[1]) {
  if (colData_MSBB_FL[i, "batch"] == "K85C014" | 
      (colData_MSBB_FL[i, "batch"] == "B18C014" & 
       !str_detect(colData_MSBB_FL[i, "Run"], "resequenced"))) {
    
    colData_MSBB_FL[i, "batch"] <- "BOTTOM"
  } 
  else if(str_detect(colData_MSBB_FL[i, "batch"], "S")) {
    colData_MSBB_FL[i, "batch"] <- "LEFT"
  }
  else {
    colData_MSBB_FL[i, "batch"] <- "RIGHT"
  }
}


# adjust ROSMAP and GSE110731 for MSBB
batch <- rep("GB", dim(colData_GSE110731)[1])
colData_GSE110731 <- cbind(colData_GSE110731, batch)

batch <- rep("RB", dim(colData_ROSMAP_FL)[1])
colData_ROSMAP_FL <- cbind(colData_ROSMAP_FL, batch)

# unify the Post-mortem interval as in minutes
colData_GSE110731[, "pmi"] <- round(as.numeric(as.character(colData_GSE110731[, "pmi"])) * 60)
colData_ROSMAP_FL[, "pmi"] <- round(as.numeric(as.character(colData_ROSMAP_FL[, "pmi"])) * 60)


colDataList_FL <- list(colData_GSE110731, colData_ROSMAP_FL, colData_MSBB_FL)
txiList_FL <- list(txi_GSE110731, txi_ROSMAP, txi_MSBB)


colTxi_FL <- Combine_ColDatas_Txis(colDataList_FL, txiList_FL, brain_region = "FL")

### without adjustment
dds_combined_FL <- ConstructDDSandVisualise(txi = colTxi_FL$txi_combined, 
                                            colData = colTxi_FL$colData_combined, 
                                            brain_region = "combined_FL", 
                                            fig_dir = "./DE_combined_test_figures/",
                                            BEC = FALSE, transform = "vst",
                                            correct_condition = FALSE,
                                            pca_intergroup = c("condition", "batch"),
                                            ggplot2_aes = aes(PC1, PC2, 
                                                              shape = condition, color = batch))


### batches + sva
dds_combined_FL_batches_sva <- ConstructDDSandVisualise(txi = colTxi_FL$txi_combined, 
                                                        colData = colTxi_FL$colData_combined, 
                                                        brain_region = "combined_FL", 
                                                        fig_dir = "./DE_combined_batches_sva_figures/",
                                                        BEC = "batches_sva", transform = "vst",
                                                        n_sv = 2,
                                                        design = ~ batch + age + pmi + condition,
                                                        correct_condition = FALSE,
                                                        pca_intergroup = c("condition", "batch"),
                                                        ggplot2_aes = aes(PC1, PC2, 
                                                                          shape = condition, 
                                                                          color = batch))


DEres_combined_FL_batches_sva <- DEanalysis(ddsM = dds_combined_FL_batches_sva$ddsM, 
                                            ddsF = dds_combined_FL_batches_sva$ddsF, 
                                            brain_region = "combined_FL", 
                                            figureSuffix = "_braak_batches_sva_vst",
                                            fig_dir = "./DE_combined_batches_sva_figures/")


### batches
dds_combined_FL_batches <- ConstructDDSandVisualise(txi = colTxi_FL$txi_combined, 
                                                    colData = colTxi_FL$colData_combined, 
                                                    brain_region = "combined_FL", 
                                                    fig_dir = "./DE_combined_batches_figures/",
                                                    BEC = "batches", transform = "vst",
                                                    design = ~ batch + age + pmi + condition,
                                                    correct_condition = FALSE,
                                                    pca_intergroup = c("condition", "batch"),
                                                    ggplot2_aes = aes(PC1, PC2, 
                                                                      shape = condition, 
                                                                      color = batch))

DEres_combined_FL_batches <- DEanalysis(ddsM = dds_combined_FL_batches$ddsM, 
                                        ddsF = dds_combined_FL_batches$ddsF, 
                                        brain_region = "combined_FL", 
                                        figureSuffix = "_braak_batches_vst",
                                        fig_dir = "./DE_combined_batches_figures/")



### sva
dds_combined_FL_sva <- ConstructDDSandVisualise(txi = colTxi_FL$txi_combined, 
                                                colData = colTxi_FL$colData_combined, 
                                                brain_region = "FL", 
                                                fig_dir = "./DE_combined_sva_figures/",
                                                BEC = "sva", transform = "vst",
                                                n_sv = 4,
                                                design = ~ condition,
                                                correct_condition = FALSE,
                                                pca_intergroup = c("condition", "batch"),
                                                ggplot2_aes = aes(PC1, PC2, 
                                                                  shape = condition, 
                                                                  color = batch))


DEres_combined_FL_sva <- DEanalysis(ddsM = dds_combined_FL_sva$ddsM, 
                                    ddsF = dds_combined_FL_sva$ddsF, 
                                    brain_region = "combined_FL", 
                                    figureSuffix = "_braak_sva_vst",
                                    fig_dir = "./DE_combined_sva_figures/")



