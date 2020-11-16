# differential expression analysis for each datasets using DESeq2
source('/binf-isilon/alab/students/vrw936/Master_Thesis/MyCode/DESeq2_functions.R')


# apply DESeq2 to each datasets to perform DE analysis and make plots ----------

setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna")


# GSE110731

txi_gse110731 <- readRDS("txi_GSE110731.rds")
colData_gse110731 <- readRDS("colData_GSE110731.rds")


dds_gse110731_sva <- ConstructDDSandVisualise(txi = txi_gse110731, 
                                              colData = colData_gse110731,
                                              brain_region = "GSE110731", transform = 'vst',
                                              design = ~ condition, 
                                              BEC = "sva", n_sv = 2,
                                              fig_dir = "./DE_sva_figures/")

DEres_gse110731_sva <- DEanalysis(ddsM = dds_gse110731_sva$ddsM, 
                                  ddsF = dds_gse110731_sva$ddsF, 
                                  brain_region = "GSE110731", 
                                  figureSuffix = "_braak_sva_vst",
                                  fig_dir = "./DE_sva_figures/")


# prepareing for RUVseq
dds_gse110731_batches <- ConstructDDSandVisualise(txi = txi_gse110731, 
                                                  colData = colData_gse110731,
                                                  brain_region = "GSE110731", 
                                                  transform = 'vst',
                                                  design = ~ age + pmi + condition,
                                                  BEC = "batches",
                                                  fig_dir = "./DE_batches_figures/")


dds_gse110731_batches_sva <- ConstructDDSandVisualise(txi = txi_gse110731, 
                                                      colData = colData_gse110731,
                                                      brain_region = "GSE110731", 
                                                      transform = 'vst',
                                                      design = ~ age + pmi + condition,
                                                      BEC = "batches_sva", n_sv = 2,
                                                      fig_dir = "./DE_batches_sva_figures/")



DEres_gse110731_batches <- DEanalysis(ddsM = dds_gse110731_batches$ddsM, 
                                      ddsF = dds_gse110731_batches$ddsF, 
                                      brain_region = "GSE110731", 
                                      figureSuffix = "_braak_batches_vst",
                                      fig_dir = "./DE_batches_figures/")

DEres_gse110731_batches_sva <- DEanalysis(ddsM = dds_gse110731_batches_sva$ddsM, 
                                          ddsF = dds_gse110731_batches_sva$ddsF, 
                                          brain_region = "GSE110731", 
                                          figureSuffix = "_braak_batches_sva_vst",
                                          fig_dir = "./DE_batches_sva_figures/")

NegCtrlGenes_gse110731 <- NegativeControlGenesFromRes(ddsM = dds_gse110731_batches_sva$ddsM,
                                                      ddsF = dds_gse110731_batches_sva$ddsF,
                                                      resM = DEres_gse110731_batches_sva$resM,
                                                      resF = DEres_gse110731_batches_sva$resF)

saveRDS(NegCtrlGenes_gse110731, file = "NegCtrlGenes_gse110731.rds")


Ws_gse110731 <- GetWsByRUVg(empiricalM = NegCtrlGenes_gse110731$empiricalM, 
                            empiricalF = NegCtrlGenes_gse110731$empiricalF,
                            ddsM = dds_gse110731_batches_sva$ddsM, 
                            resM = DEres_gse110731_batches_sva$resM, 
                            ddsF = dds_gse110731_batches_sva$ddsF, 
                            resF = DEres_gse110731_batches_sva$resF)


dds_gse110731_ruvseq <- ConstructDDSandVisualise(txi = txi_gse110731, 
                                                 colData = colData_gse110731,
                                                 brain_region = "GSE110731", 
                                                 transform = 'vst',
                                                 design = ~ condition, 
                                                 BEC = "RUVseq", 
                                                 Ws = Ws_gse110731,
                                                 fig_dir = "./DE_ruvseq_figures/")

DEres_gse110731_ruvseq <- DEanalysis(ddsM = dds_gse110731_ruvseq$ddsM, 
                                     ddsF = dds_gse110731_ruvseq$ddsF,
                                     brain_region = "GSE110731", 
                                     figureSuffix = "_braak_ruvseq_vst",
                                     fig_dir = "./DE_ruvseq_figures/")

# GSE53697 (sex of all samples to be predicted)

txi_gse53697 <- readRDS("txi_GSE53697.rds")
colData_gse53697 <- readRDS("colData_GSE53697.rds")

dds_gse53697_sva <- ConstructDDSandVisualise(txi = txi_gse53697, colData = colData_gse53697,
                                              brain_region = "GSE53697", transform = 'vst',
                                              design = ~ condition, BEC = "sva")


DEres_gse53697_sva <- DEanalysis(ddsM = dds_gse53697_sva$ddsM, ddsF = dds_gse53697_sva$ddsF, 
                                  brain_region = "GSE53697", figureSuffix = "_braak_sva_vst")


# GSE104704 (sex of all samples to be predicted)
txi_gse104704 <- readRDS("txi_GSE104704.rds")
colData_gse104704 <- readRDS("colData_GSE104704.rds")


dds_gse104704_sva <- ConstructDDSandVisualise(txi = txi_gse104704, colData = colData_gse104704,
                                             brain_region = "GSE104704", transform = 'vst',
                                             design = ~ condition, BEC = "sva")

DEres_gse104704_sva <- DEanalysis(ddsM = dds_gse104704_sva$ddsM, ddsF = dds_gse104704_sva$ddsF, 
                                 brain_region = "GSE104704", figureSuffix = "_braak_sva_vst")


# GSE125583 
txi_gse125583 <- readRDS("txi_GSE125583.rds")
colData_gse125583 <- readRDS("colData_GSE125583.rds")

dds_gse125583_sva <- ConstructDDSandVisualise(txi = txi_gse125583, colData = colData_gse125583,
                                              brain_region = "GSE125583", transform = 'vst',
                                              design = ~ condition, BEC = "sva")

DEres_gse125583_sva <- DEanalysis(ddsM = dds_gse125583_sva$ddsM, ddsF = dds_gse125583_sva$ddsF, 
                                  brain_region = "GSE125583", figureSuffix = "_braak_sva_vst")



# GSE95587
txi_gse95587 <- readRDS("txi_GSE95587.rds")
colData_gse95587 <- readRDS("colData_GSE95587.rds")

dds_gse95587_sva <- ConstructDDSandVisualise(txi = txi_gse95587, colData = colData_gse95587,
                                              brain_region = "GSE95587", transform = 'vst',
                                              design = ~ condition, BEC = "sva")

DEres_gse95587_sva <- DEanalysis(ddsM = dds_gse95587_sva$ddsM, ddsF = dds_gse95587_sva$ddsF, 
                                  brain_region = "GSE95587", figureSuffix = "_braak_sva_vst")


# ROSMAP (FL)
txi_ROSMAP_FL <- readRDS("txi_ROSMAP.rds")
colData_ROSMAP_FL <- readRDS("colData_ROSMAP_FL.rds")


dds_ROSMAP_FL_sva <- ConstructDDSandVisualise(txi = txi_ROSMAP_FL, 
                                              colData = colData_ROSMAP_FL,
                                              brain_region = "ROSMAP_FL", transform = 'vst',
                                              design = ~ condition, 
                                              BEC = "sva", n_sv = 2)

DEres_ROSMAP_FL_sva <- DEanalysis(ddsM = dds_ROSMAP_FL_sva$ddsM, 
                                  ddsF = dds_ROSMAP_FL_sva$ddsF, 
                                  brain_region = "ROSMAP_FL", 
                                  figureSuffix = "_braak_sva_vst")


dds_ROSMAP_FL_batches_sva <- ConstructDDSandVisualise(txi = txi_ROSMAP_FL, 
                                                      colData = colData_ROSMAP_FL,
                                                      brain_region = "ROSMAP_FL", 
                                                      transform = 'vst',
                                                      design = ~ age + pmi + condition,
                                                      BEC = "batches_sva",
                                                      n_sv = 2,
                                                      fig_dir = "./DE_batches_sva_figures/")

DEres_ROSMAP_FL_batches_sva <- DEanalysis(ddsM = dds_ROSMAP_FL_batches_sva$ddsM, 
                                          ddsF = dds_ROSMAP_FL_batches_sva$ddsF, 
                                          brain_region = "ROSMAP_FL", 
                                          figureSuffix = "_braak_batches_sva_vst",
                                          fig_dir = "./DE_batches_sva_figures/")

NegCtrlGenes_ROSMAP_FL <- NegativeControlGenesFromRes(ddsM = dds_ROSMAP_FL_batches_sva$ddsM,
                                                      ddsF = dds_ROSMAP_FL_batches_sva$ddsF,
                                                      resM = DEres_ROSMAP_FL_batches_sva$resM,
                                                      resF = DEres_ROSMAP_FL_batches_sva$resF)

saveRDS(NegCtrlGenes_ROSMAP_FL, file = "NegCtrlGenes_ROSMAP_FL.rds")


Ws_ROSMAP_FL <- GetWsByRUVg(empiricalM = NegCtrlGenes_ROSMAP_FL$empiricalM, 
                            empiricalF = NegCtrlGenes_ROSMAP_FL$empiricalF,
                            ddsM = dds_ROSMAP_FL_batches_sva$ddsM, 
                            resM = DEres_ROSMAP_FL_batches_sva$resM, 
                            ddsF = dds_ROSMAP_FL_batches_sva$ddsF, 
                            resF = DEres_ROSMAP_FL_batches_sva$resF)


dds_ROSMAP_FL_ruvseq <- ConstructDDSandVisualise(txi = txi_ROSMAP_FL, 
                                                 colData = colData_ROSMAP_FL,
                                                 brain_region = "ROSMAP_FL", 
                                                 figureSuffix = "_braak_ruvseq_vst",
                                                 fig_dir = "./DE_ruvseq_figures/")




# MSBB FL
txi_MSBB <- readRDS("txi_MSBB.rds")
colData_MSBB_FL <- readRDS("colData_MSBB_FL.rds")


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


dds_MSBB_FL <- ConstructDDSandVisualise(txi = txi_MSBB, colData = colData_MSBB_FL,
                                            brain_region = "MSBB_FL", transform = 'vst',
                                            design = ~ condition, BEC = FALSE, 
                                            pca_intergroup = c("batch", "condition"),
                                            ggplot2_aes = aes(PC1, PC2, color = batch),
                                            fig_dir = "./DE_test_figures/")

dds_MSBB_FL_sva <- ConstructDDSandVisualise(txi = txi_MSBB, colData = colData_MSBB_FL,
                                            brain_region = "MSBB_FL", transform = 'vst',
                                            design = ~ condition, BEC = "sva", 
                                            n_sv = 2,
                                            pca_intergroup = c("batch", "condition"),
                                            ggplot2_aes = aes(PC1, PC2, color = batch),
                                            fig_dir = "./DE_sva_figures/")

DEres_MSBB_FL_sva <- DEanalysis(ddsM = dds_MSBB_FL_sva$ddsM, ddsF = dds_MSBB_FL_sva$ddsF, 
                                brain_region = "MSBB_FL", figureSuffix = "_braak_sva_vst")

dds_MSBB_FL_batches <- ConstructDDSandVisualise(txi = txi_MSBB, 
                                                colData = colData_MSBB_FL,
                                                brain_region = "MSBB_FL", 
                                                transform = 'vst',
                                                design = ~ batch + age + pmi + condition,
                                                BEC = "batches",
                                                pca_intergroup = c("batch", "condition"),
                                                ggplot2_aes = aes(PC1, PC2, color = batch),
                                                fig_dir = "./DE_batches_figures/")


DEres_MSBB_FL_batches <- DEanalysis(ddsM = dds_MSBB_FL_batches$ddsM, 
                                    ddsF = dds_MSBB_FL_batches$ddsF, 
                                    brain_region = "MSBB_FL", 
                                    figureSuffix = "_braak_batches_vst",
                                    fig_dir = "./DE_batches_figures/")



dds_MSBB_FL_batches_sva <- ConstructDDSandVisualise(txi = txi_MSBB, 
                                                    colData = colData_MSBB_FL,
                                                    brain_region = "MSBB_FL", 
                                                    transform = 'vst',
                                                    design = ~ batch + age + pmi + condition,
                                                    BEC = "batches_sva", n_sv = 2,
                                                    pca_intergroup = c("batch", "condition"),
                                                    ggplot2_aes = aes(PC1, PC2, shape = batch, 
                                                                      color = condition),
                                                    fig_dir = "./DE_batches_sva_figures/")

DEres_MSBB_FL_batches_sva <- DEanalysis(ddsM = dds_MSBB_FL_batches_sva$ddsM, 
                                        ddsF = dds_MSBB_FL_batches_sva$ddsF, 
                                        brain_region = "MSBB_FL", 
                                        figureSuffix = "_braak_batches_sva_vst",
                                        fig_dir = "./DE_batches_sva_figures/")

NegCtrlGenes_MSBB_FL <- NegativeControlGenesFromRes(ddsM = dds_MSBB_FL_batches_sva$ddsM,
                                                    ddsF = dds_MSBB_FL_batches_sva$ddsF,
                                                    resM = DEres_MSBB_FL_batches_sva$resM,
                                                    resF = DEres_MSBB_FL_batches_sva$resF)

saveRDS(NegCtrlGenes_MSBB_FL, file = "NegCtrlGenes_MSBB_FL.rds")


Ws_MSBB_FL <- GetWsByRUVg(empiricalM = NegCtrlGenes_MSBB_FL$empiricalM, 
                          empiricalF = NegCtrlGenes_MSBB_FL$empiricalF,
                          ddsM = dds_MSBB_FL_batches_sva$ddsM, 
                          resM = DEres_MSBB_FL_batches_sva$resM, 
                          ddsF = dds_MSBB_FL_batches_sva$ddsF, 
                          resF = DEres_MSBB_FL_batches_sva$resF)


dds_MSBB_FL_ruvseq <- ConstructDDSandVisualise(txi = txi_MSBB_FL, 
                                               colData = colData_MSBB_FL,
                                               brain_region = "MSBB_FL", 
                                               figureSuffix = "_braak_ruvseq_vst",
                                               fig_dir = "./DE_ruvseq_figures/")





# MSBB TL
colData_MSBB_TL <- readRDS("colData_MSBB_TL.rds")

dds_MSBB_TL_sva <- ConstructDDSandVisualise(txi = txi_MSBB, colData = colData_MSBB_TL,
                                            brain_region = "MSBB_TL", transform = 'vst',
                                            design = ~ condition, BEC = "sva")

DEres_MSBB_TL_sva <- DEanalysis(ddsM = dds_MSBB_TL_sva$ddsM, ddsF = dds_MSBB_TL_sva$ddsF, 
                                brain_region = "MSBB_TL", figureSuffix = "_braak_sva_vst")


# MayoRNAseq CB

txi_MayoRNAseq_CBE <- readRDS("txi_MayoRNAseq_CB.rds")
colData_MayoRNAseq_CBE <- readRDS("colData_MayoRNAseq_CBE.rds")


dds_MayoRNAseq_CBE_sva <- ConstructDDSandVisualise(txi = txi_MayoRNAseq_CBE, 
                                                   colData = colData_MayoRNAseq_CBE,
                                                   brain_region = "MayoRNAseq_CBE", 
                                                   transform = 'vst',
                                                   design = ~ condition, BEC = "sva")

DEres_MayoRNAseq_CBE_sva <- DEanalysis(ddsM = dds_MayoRNAseq_CBE_sva$ddsM, 
                                       ddsF = dds_MayoRNAseq_CBE_sva$ddsF,
                                       brain_region = "MayoRNAseq_CBE", 
                                       figureSuffix = "_braak_sva_vst")


# MayoRNAseq TL

txi_MayoRNAseq_TCX <- readRDS("txi_MayoRNAseq_TCX.rds")
colData_MayoRNAseq_TCX <- readRDS("colData_MayoRNAseq_TCX.rds")

dds_MayoRNAseq_TCX_sva <- ConstructDDSandVisualise(txi = txi_MayoRNAseq_TCX, 
                                                   colData = colData_MayoRNAseq_TCX,
                                                   brain_region = "MayoRNAseq_TCX", 
                                                   transform = 'vst',
                                                   design = ~ condition, BEC = "sva")

DEres_MayoRNAseq_TCX_sva <- DEanalysis(ddsM = dds_MayoRNAseq_TCX_sva$ddsM, 
                                       ddsF = dds_MayoRNAseq_TCX_sva$ddsF,
                                       brain_region = "MayoRNAseq_TCX", 
                                       figureSuffix = "_braak_sva_vst")




# bar plot and Venn plot for each brain region ----------------------------


BarplotVennplot <- function(res_dir = "DE_batches_sva_figures",
                            fig_dir = "./Bar_Venn/", brain_region = "FL",
                            RegionKeys = c("ROSMAP_FL", "GSE110731", "MSBB_FL"),
                            colours = c("yellow", "blue", "green"),
                            alpha = 0.05){
  # make barplots and Venn plots showing the the number of DEGs
  # comparing different datasets and sexes, for a specified brain region
  
  dir.create(fig_dir)
  
  ##### a bar plot for the brain region, including all datasets
  
  # make a 'gene - sex - dataset' tibble combining all datasets and sexes for the brain region
  # preparing for bar plot
  
  # load DE results
  for (DEres in list.files(path = res_dir, pattern = "DEres")) {
    assign(gsub("(^.+)\\.rds", "\\1", DEres), readRDS(file.path(res_dir, DEres)))
  }
  
  res_combined <- tibble(gene = character(), sex = character(), dataset = character())
  
  for (res in ls(pattern = "DEres.+batches_sva")){  # M and F are separated
    sex_info <- gsub("DEres_.+_([F|M])_.+", "\\1", res)
    dataset_info <- gsub("DEres_(.+)_[F|M]_.+", "\\1", res)
    
    if(dataset_info %in% RegionKeys){
      
      assign(res, dplyr::filter(as_tibble(get(res), rownames = "gene"), 
                                padj < alpha))
      assign(res, dplyr::select(get(res), gene))
      sex_dataset <- tibble(sex = rep(sex_info, dim(get(res))[1]), 
                            dataset = rep(dataset_info, dim(get(res))[1]))
      assign(res, bind_cols(get(res), sex_dataset))
      
      res_combined <- bind_rows(res_combined, get(res))
    }
  }
  
  
  ggplot(data = res_combined, mapping = aes(x = dataset, fill = sex)) +
    geom_bar(position = "dodge") +
    ggsave(paste0(fig_dir, "barplot_", brain_region, ".png"))
  
  
  ##### 2 Venn plot for Male and Female respectively, with each across all datasets
  
  # load the data again since the original ones have been changed 
  for (DEres in list.files(path = res_dir, pattern = "DEres")) {
    assign(gsub("(^.+)\\.rds", "\\1", DEres), readRDS(file.path(res_dir, DEres)))
  }
  
  M_list <- list()
  M_names <- c()
  F_list <- list()
  F_names <- c()
  
  for (res in ls(pattern = "DEres.+batches_sva")){  # M and F are separated
    sex_info <- gsub("DEres_.+_([F|M])_.+", "\\1", res)
    dataset_info <- gsub("DEres_(.+)_[F|M]_.+", "\\1", res)
    
    if(dataset_info %in% RegionKeys){
      assign(res, filter(as_tibble(get(res), rownames = "gene"), 
                         padj < alpha))
      
      if(sex_info == "M"){
        M_names <- c(M_names, dataset_info)
        
        tmp_list_M <- list(get(res)$gene)
        M_list <- c(M_list, tmp_list_M)
      }
      else if(sex_info == "F"){
        F_names <- c(F_names, dataset_info)
        tmp_list_F <- list(get(res)$gene)
        F_list <- c(F_list, tmp_list_F)
      }
    }
  }
  
  names(M_list) <- M_names
  names(F_list) <- F_names
  
  venn.diagram(x = M_list, 
               filename = file.path(fig_dir, paste0("Venn_", brain_region, "_M.tiff")),
               fill = colours, 
               cat.cex = 1, cex = 3, width = 4000)
  
  venn.diagram(x = F_list, 
               filename = file.path(fig_dir, paste0("Venn_", brain_region, "_F.tiff")), 
               fill = colours, 
               cat.cex = 1, cex = 3, width = 4000)
  
  
  ##### one Venn plot for each dataset, across male and female
  
  # load the data again since the original ones have been changed
  for (DEres in list.files(path = res_dir, pattern = "DEres")) {
    assign(gsub("(^.+)\\.rds", "\\1", DEres), readRDS(file.path(res_dir, DEres)))
  }   
  
  for(dataset in RegionKeys){
    VennX = list()
    Names_sex = c()
    
    for (res in ls(pattern = "DEres.+batches_sva")){  # M and F are separated
      sex_info <- gsub("DEres_.+_([F|M])_.+", "\\1", res)
      dataset_info <- gsub("DEres_(.+)_[F|M]_.+", "\\1", res)
      
      if(dataset_info == dataset){
        assign(res, filter(as_tibble(get(res), rownames = "gene"), 
                           padj < alpha))
        
        if(sex_info == "M"){
          tmp_list_X <- list(get(res)$gene)
          VennX <- c(VennX, tmp_list_X)
          Names_sex <- c(Names_sex, "M")
        }
        else if(sex_info == "F"){
          tmp_list_X <- list(get(res)$gene)
          VennX <- c(VennX, tmp_list_X)
          Names_sex <- c(Names_sex, "F")
        }
      }
    }
    
    names(VennX) <- Names_sex
    
    venn.diagram(x = VennX, 
                 filename = file.path(fig_dir, paste0("Venn_", dataset, "_.tiff")),
                 fill = c("cornflowerblue", "darkorchid1"), 
                 cat.cex = 1, cex = 3, width = 4000)
  }
  
  system(paste("rm", file.path(fig_dir, "*.log")))
}



setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/")

### frontal lobe (FL)

BarplotVennplot(res_dir = "DE_batches_sva_figures",
                fig_dir = "./Bar_Venn/", brain_region = "FL",
                RegionKeys = c("ROSMAP_FL", "GSE110731", "MSBB_FL"))

DEres_combined_FL_F_braak_batches_sva_vst <- readRDS("DE_combined_batches_sva_figures/DEres_combined_FL_F_braak_batches_sva_vst.rds")
DEres_combined_FL_M_braak_batches_sva_vst <- readRDS("DE_combined_batches_sva_figures/DEres_combined_FL_M_braak_batches_sva_vst.rds")

BarplotVennplot(res_dir = "DE_batches_sva_figures",
                fig_dir = "./Bar_Venn/", brain_region = "FL",
                RegionKeys = c("ROSMAP_FL", "GSE110731", "MSBB_FL", "combined_FL"), 
                colours = c("yellow", "blue", "green", "red"))



