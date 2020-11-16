# PathwayTopoAnalysisVisual.R
# analyse and visualise the topology of pathway
library(graphite)
library(pathview)
library(tidyverse)
library(stringr)
library(clusterProfiler)


# pathview for KEGG pathways ----------------------------------------------

PathwayVisualise <- function(res, fig_dir = "Pathway_figures/", brain_region,
                             alpha = 0.05, lfc_threshold = 0.1){
  # a plot showing the topological structure and log2FC value of member DEGs
  # for each enriched pathway
  # pathways in KEGG are visualised by pathview
  # other pathways maybe visualised by ToPASeq and Graphite
  
  
  dir.create(fig_dir)
  setwd(fig_dir)
  
  
  res <- filter(as_tibble(res, rownames = "gene"), 
                padj < alpha, abs(log2FoldChange) > lfc_threshold)
  symbol_entrezid <- bitr(res$gene, fromType="SYMBOL", toType="ENTREZID",
                          OrgDb="org.Hs.eg.db")
  colnames(symbol_entrezid) <- c("gene", "entrezid")
  res_entrezid <- inner_join(res, symbol_entrezid, by = "gene")
  geneList <- res_entrezid$log2FoldChange
  names(geneList) <- res_entrezid$entrezid
  
  geneList <- sort.default(geneList, decreasing = TRUE)  # must use sort.default(), not sort()
  gene <- names(geneList)  # set no threshold on fold change
  
  # Mitophagy - Animal
  pathview(gene.data = geneList, 
           limit = list(gene = 1.1, cpd = 1.1),  # for log2FC
           bins = list(gene = 22, cpd = 22),  # also for colour legend
           gene.idtype = "entrez", 
           # kegg.dir = fig_dir,  # provide the xml file by myself in this dir
           pathway.id = "04137", species = "hsa",
           same.layer = FALSE, kegg.native = FALSE, 
           out.suffix = brain_region)
  
  system(paste0("mv ", "hsa04137.", brain_region, ".pdf ", 
                "kegg_mitophagy_", brain_region, ".pdf"))
  
  # Alzheimer's disease
  pathview(gene.data = geneList,
           limit = list(gene = 1.1, cpd = 1.1),  # for log2FC
           bins = list(gene = 22, cpd = 22),  # also for colour legend
           gene.idtype = "entrez", 
           # kegg.dir = fig_dir,  # provide the xml file by myself in this dir
           pathway.id = "05010", species = "hsa",
           same.layer = FALSE, kegg.native = FALSE, pdf.size = c(10, 10),
           out.suffix = brain_region)
  
  system(paste0("mv ", "hsa05010.", brain_region, ".pdf ", 
                "kegg_alzheimer_", brain_region, ".pdf"))
  
  # Lysosome
  pathview(gene.data = geneList,
           limit = list(gene = 1.1, cpd = 1.1),  # for log2FC
           bins = list(gene = 22, cpd = 22),  # also for colour legend
           gene.idtype = "entrez", 
           # kegg.dir = fig_dir,  # provide the xml file by myself in this dir
           pathway.id = "04142", species = "hsa",
           same.layer = TRUE, kegg.native = TRUE, 
           out.suffix = brain_region)
  
  system(paste0("mv ", "hsa04142.", brain_region, ".png ", 
                "kegg_lysosome_", brain_region, ".png"))
  
  # Base excision repair
  pathview(gene.data = geneList,
           limit = list(gene = 1.1, cpd = 1.1),  # for log2FC
           bins = list(gene = 22, cpd = 22),  # also for colour legend
           gene.idtype = "entrez", 
           # kegg.dir = fig_dir,  # provide the xml file by myself in this dir
           pathway.id = "03410", species = "hsa",
           same.layer = TRUE, kegg.native = TRUE, 
           out.suffix = brain_region)
  
  system(paste0("mv ", "hsa03410.", brain_region, ".png ", 
                "kegg_base_excision_repair_", brain_region, ".png"))
  
  
  # system(paste("mv", "hsa*", fig_dir))
  system(paste("rm", "hsa*"))
  setwd("..")
  
}




setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/")

dds_combined_FL_F_batches_sva_vst <- readRDS("DE_combined_batches_sva_figures/dds_combined_FL_F_batches_sva_vst.rds")
dds_combined_FL_M_batches_sva_vst <- readRDS("DE_combined_batches_sva_figures/dds_combined_FL_M_batches_sva_vst.rds")
DEres_combined_FL_F_braak_batches_sva_vst <- readRDS("DE_combined_batches_sva_figures/DEres_combined_FL_F_braak_batches_sva_vst.rds")
DEres_combined_FL_M_braak_batches_sva_vst <- readRDS("DE_combined_batches_sva_figures/DEres_combined_FL_M_braak_batches_sva_vst.rds")

PathwayVisualise(res = DEres_combined_FL_M_braak_batches_sva_vst, 
                 brain_region = "combined_FL_M", alpha = 0.05)

PathwayVisualise(res = DEres_combined_FL_F_braak_batches_sva_vst, 
                 brain_region = "combined_FL_F", alpha = 0.05)


# Graphite for other pathway with known structure -------------------------


Graphite_Cytoscape <- function(fig_dir = "Pathway_figures/", 
                               key_word, DE_result, alpha = 0.05,
                               lfc_threshold = 0.1, brain_region){
  
  DE_result <- filter(as_tibble(DE_result, rownames = "gene"), 
                padj < alpha, abs(log2FoldChange) > lfc_threshold)
  
  pathDBs <- pathwayDatabases()
  pathDBs <- filter(pathDBs, species == "hsapiens")
  
  for (j in 1:dim(pathDBs)[1]) {
    pathDB <- pathways(pathDBs$species[j], pathDBs$database[j])
    path_names_tmp <- names(pathDB)[str_detect(names(pathDB), 
                                               regex(key_word, ignore_case = TRUE))]
    print(as.character(pathDBs$database[j]))
    print(path_names_tmp)
    
    
    if (length(path_names_tmp) != 0) {
      for (path_name in path_names_tmp) {
        fig_name <- gsub("[/  ]", "_", tolower(path_name))
        fig_name <- paste0(as.character(pathDBs$database[j]), "_", fig_name)
        if(fig_name == "panther_oxidative_stress_response"){next}  # due to wrong id
        if(fig_name == "panther_alzheimer_disease-amyloid_secretase_pathway"){next}
        if(fig_name == "panther_alzheimer_disease-presenilin_pathway"){next}
        # NOTE: these pathways (incl. others in panther) have UNIPROT ID of many species, 
        # including rat, mouse, zebra fish, etc. which cannot be converted using bitr(OrgDb = "org.Hs.eg.db")
        # so they cannot be converted to gene symbols
        
        pathway <- pathDB[[path_name]]
        cytoscape_table <- edges(pathway)
        id_type <- as.character(cytoscape_table$src_type[1])
        cytoscape_table$src_type <- 0
        cytoscape_table$dest_type <- 0
         
        for(i in 1:dim(cytoscape_table)[1]){
          src_tmp <- bitr(cytoscape_table$src[i], 
                          fromType = id_type, 
                          toType = "SYMBOL",
                          OrgDb = "org.Hs.eg.db")[, 2]
          if(!is.na(src_tmp)){cytoscape_table$src[i] <- src_tmp}else{print("unmappable")}
          
          dest_tmp <- bitr(cytoscape_table$dest[i], 
                           fromType = id_type, 
                           toType = "SYMBOL",
                           OrgDb = "org.Hs.eg.db")[, 2]
          if(!is.na(dest_tmp)){cytoscape_table$dest[i] <- dest_tmp}else{print("unmappable")}
          
          
          if(cytoscape_table$src[i] %in% DE_result$gene){
            cytoscape_table$src_type[i] <- DE_result[DE_result$gene == 
                                                       cytoscape_table$src[i], ]$log2FoldChange
          }
          
          if(cytoscape_table$dest[i] %in% DE_result$gene){
            cytoscape_table$dest_type[i] <- DE_result[DE_result$gene == 
                                                        cytoscape_table$dest[i], ]$log2FoldChange
          }
          
        }
        
        write_tsv(x = cytoscape_table, path = file.path(fig_dir,
                                                        paste0(fig_name, "_",
                                                               brain_region, ".tsv")))
      }
    }
  }
}



# mitophagy
Graphite_Cytoscape(fig_dir = "Pathway_figures/",
                   key_word = "mitophag", 
                   DE_result = DEres_combined_FL_M_braak_batches_sva_vst, 
                   alpha = 0.05,
                   lfc_threshold = 0.1, 
                   brain_region = "combined_FL_M")


Graphite_Cytoscape(fig_dir = "Pathway_figures/",
                   key_word = "mitophag", 
                   DE_result = DEres_combined_FL_F_braak_batches_sva_vst, 
                   alpha = 0.05,
                   lfc_threshold = 0.1, 
                   brain_region = "combined_FL_F")



# lysosome
Graphite_Cytoscape(fig_dir = "Pathway_figures/",
                   key_word = "lysoso", 
                   DE_result = DEres_combined_FL_M_braak_batches_sva_vst, 
                   alpha = 0.05,
                   lfc_threshold = 0.1, 
                   brain_region = "combined_FL_M")


Graphite_Cytoscape(fig_dir = "Pathway_figures/",
                   key_word = "lysoso", 
                   DE_result = DEres_combined_FL_F_braak_batches_sva_vst, 
                   alpha = 0.05,
                   lfc_threshold = 0.1, 
                   brain_region = "combined_FL_F")



# oxidative stress response 

Graphite_Cytoscape(fig_dir = "Pathway_figures/",
                   key_word = "oxidative stress", 
                   DE_result = DEres_combined_FL_M_braak_batches_sva_vst, 
                   alpha = 0.05,
                   lfc_threshold = 0.1, 
                   brain_region = "combined_FL_M")


Graphite_Cytoscape(fig_dir = "Pathway_figures/",
                   key_word = "oxidative stress", 
                   DE_result = DEres_combined_FL_F_braak_batches_sva_vst, 
                   alpha = 0.05,
                   lfc_threshold = 0.1, 
                   brain_region = "combined_FL_F")



# Alzheimer
Graphite_Cytoscape(fig_dir = "Pathway_figures/",
                   key_word = "Alzheimer", 
                   DE_result = DEres_combined_FL_M_braak_batches_sva_vst, 
                   alpha = 0.05,
                   lfc_threshold = 0.1, 
                   brain_region = "combined_FL_M")


Graphite_Cytoscape(fig_dir = "Pathway_figures/",
                   key_word = "Alzheimer", 
                   DE_result = DEres_combined_FL_F_braak_batches_sva_vst, 
                   alpha = 0.05,
                   lfc_threshold = 0.1, 
                   brain_region = "combined_FL_F")



##### bar plot showing DEG's lfc 

for(tsv in list.files(path = "Pathway_figures", pattern = "tsv")){
  fig_name = gsub("(^.+)\\.tsv", "\\1", tsv)
  
  tsv_tmp <- read_tsv(file.path("Pathway_figures", tsv))
  tsv_tmp <- bind_rows(dplyr::select(tsv_tmp, symbol = src, log2FoldChange = src_type), 
                       dplyr::select(tsv_tmp, symbol = dest, log2FoldChange = dest_type))
  
  tsv_tmp <- unique.data.frame(tsv_tmp)
  tsv_tmp <- dplyr::filter(tsv_tmp, log2FoldChange != 0)
  
  ggplot(data = tsv_tmp, mapping = aes(x = reorder(symbol, log2FoldChange), 
                                       y = log2FoldChange)) +
    geom_bar(stat = "identity") +
    labs(x = "Gene") +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.7)) +
    ggsave(filename = paste0(fig_name, "_lfc.png"), path = "Pathway_figures")
  
}
