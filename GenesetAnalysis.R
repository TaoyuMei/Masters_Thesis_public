# geneset analysis

library(tidyverse)
library(pathview)
library(clusterProfiler)
library(msigdbr)
library(GSVA)
library(ReactomePA)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(xtable)
library(ggpubr)


# prepare data for gene set analysis and visualisation --------------------------------------


##### wikiPathways

wp2gene <- read.gmt("/binf-isilon/alab/students/vrw936/scratch/reference_ome_and_anno/wikipathways-20200510-gmt-Homo_sapiens.gmt")
# http://data.wikipathways.org/current/gmt/wikipathways-20200510-gmt-Homo_sapiens.gmt
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

##### MSigDB
msigdbr <- msigdbr(species = "Homo sapiens")
ms2gene <- dplyr::select(msigdbr, gs_name, entrez_gene)


##### EntrezID - Gene symbol conversion

setwd("/binf-isilon/alab/students/vrw936/scratch/reference_ome_and_anno/")
gencode.v33.metadata.HGNC <- read_tsv("gencode.v33.metadata.HGNC",
                                      col_names = FALSE)
gencode.v33.metadata.EntrezGene <- read_tsv("gencode.v33.metadata.EntrezGene",
                                            col_names = FALSE)
symbol_entrezid <- inner_join(gencode.v33.metadata.HGNC, gencode.v33.metadata.EntrezGene, by = "X1")

symbol_entrezid <- unique.data.frame(dplyr::select(symbol_entrezid, symbol = X2.x, entrezid = X2.y))


# functionalise all gene set analysis (diff. enrichment methods & databases) -------------------


DetermineBGgenes <- function(dds, symbol_entrezid){
  # use the combined dds of a brain region to determain the background genes 
  # that are actually expressed in this brain region
  
  dds <- estimateSizeFactors(dds)  # prepare for normalisation
  dat  <- counts(dds, normalized = TRUE)
  idx  <- rowMeans(dat) > 1  
  dat  <- dat[idx, ]
  bg <- tibble(SYMBOL = rownames(dat))
  # bg <- symbol_entrezid[symbol_entrezid$symbol %in% bg, ]$entrezid
  
  symbol_entrezid <- bitr(bg$SYMBOL, fromType="SYMBOL", toType="ENTREZID",
                          OrgDb="org.Hs.eg.db")
  bg_entrezid <- inner_join(bg, symbol_entrezid, by = "SYMBOL")
  
  # bg <- bitr(bg, fromType="SYMBOL", toType="ENTREZID",
  #            OrgDb="org.Hs.eg.db")[, 2]
  
  return(bg_entrezid$ENTREZID)
}


GeneSetAnalysis <- function(res, bg_gene, wpid2gene, wpid2name, ms2gene,
                            fig_dir = "GeneSet_figures", brain_region,
                            alpha = 0.05, lfc_threshold = 0.1){
  # input 1: an res object from DESeq2, containing a DEG list(incl. up and down regulated genes)
  # and relevant values (log2FC) for each gene
  # input 2: background genes
  # input 3 & 4: for WikiPathways analysis
  # input 5: for MSigDb
  # output: a "gene set name - member DEG" table for each gene set database
  
  system("date")
  dir.create(fig_dir)
  
  ##### prepare 'geneList' and 'gene'
  
  
  res <- filter(as_tibble(res, rownames = "gene"), 
                padj < alpha, abs(log2FoldChange) > lfc_threshold)  # default alpha = 0.05
  symbol_entrezid <- bitr(res$gene, fromType="SYMBOL", toType="ENTREZID",
                          OrgDb="org.Hs.eg.db")
  colnames(symbol_entrezid) <- c("gene", "entrezid")
  res_entrezid <- inner_join(res, symbol_entrezid, by = "gene")
  geneList <- 2 ** res_entrezid$log2FoldChange  # log2FC to FC
  names(geneList) <- res_entrezid$entrezid
  
  # split up- and down- regulated genes
  
  geneList_up <- sort.default(geneList[geneList > 1], decreasing = TRUE)
  geneList_down <- sort.default(geneList[geneList < 1], decreasing = TRUE)
  
  gene_up <- names(geneList_up)
  gene_down <- names(geneList_down)
  
  geneList <- sort.default(geneList, decreasing = TRUE)
  gene <- names(geneList)  
 
  
  
  ##### overrepresentation analysis
  
  # if only need to map DEGs to gene sets regardless of significance, 
  # just set pvalueCutoff = 1, qvalueCutoff = 1
  
  # up 
  
  GO_up <- enrichGO(gene = gene_up, OrgDb = org.Hs.eg.db, universe = bg_gene,
                    pvalueCutoff = 1, qvalueCutoff = 1)
  KEGGpathway_up <- enrichKEGG(gene = gene_up, keyType = "kegg", universe = bg_gene,
                               pvalueCutoff = 1, qvalueCutoff = 1)
  KEGGmodule_up <- enrichMKEGG(gene = gene_up, keyType = "kegg", universe = bg_gene,
                               pvalueCutoff = 1, qvalueCutoff = 1)
  
  wikiPathway_up <- enricher(gene_up, TERM2GENE = wpid2gene, TERM2NAME = wpid2name,
                             pvalueCutoff = 1, qvalueCutoff = 1,
                             universe = bg_gene)
  
  msig_up <- enricher(gene_up, TERM2GENE = ms2gene, universe = bg_gene,
                      pvalueCutoff = 1, qvalueCutoff = 1)
  
  reactome_up <- ReactomePA::enrichPathway(gene = gene_up, universe = bg_gene,
                                           pvalueCutoff = 1, qvalueCutoff = 1)
  
  # down
  
  GO_down <- enrichGO(gene = gene_down, OrgDb = org.Hs.eg.db, universe = bg_gene,
                      pvalueCutoff = 1, qvalueCutoff = 1)
  KEGGpathway_down <- enrichKEGG(gene = gene_down, keyType = "kegg", universe = bg_gene, 
                                 pvalueCutoff = 1, qvalueCutoff = 1)
  KEGGmodule_down <- enrichMKEGG(gene = gene_down, keyType = "kegg", universe = bg_gene,
                                 pvalueCutoff = 1, qvalueCutoff = 1)
  
  wikiPathway_down <- enricher(gene_down, TERM2GENE = wpid2gene, TERM2NAME = wpid2name,
                               pvalueCutoff = 1, qvalueCutoff = 1,
                               universe = bg_gene)
  
  msig_down <- enricher(gene_down, TERM2GENE = ms2gene, universe = bg_gene,
                        pvalueCutoff = 1, qvalueCutoff = 1)
  
  reactome_down <- ReactomePA::enrichPathway(gene = gene_down, universe = bg_gene,
                                             pvalueCutoff = 1, qvalueCutoff = 1)
  
  
  geneset_list_up <- list(GO_up = GO_up, 
                          KEGGpathway_up = KEGGpathway_up, 
                          KEGGmodule_up = KEGGmodule_up, 
                          wikiPathway_up = wikiPathway_up, 
                          msig_up = msig_up, 
                          reactome_up = reactome_up)
  
  saveRDS(geneset_list_up, file = file.path(fig_dir, paste0(brain_region, 
                                                         "_geneset_list_up.rds")))
  
  geneset_list_down <- list(GO_down = GO_down, 
                            KEGGpathway_down = KEGGpathway_down, 
                            KEGGmodule_down = KEGGmodule_down, 
                            wikiPathway_down = wikiPathway_down, 
                            msig_down = msig_down, 
                            reactome_down = reactome_down)
  
  saveRDS(geneset_list_down, file = file.path(fig_dir, paste0(brain_region, 
                                                         "_geneset_list_down.rds")))
  
  system("date")
  
}



# visualise the gene expression changes -----------------------------------

VisualiseGeneChange_keywords <- function(key_word, fig_name, symbol_entrezid, 
                                         fig_dir = "GeneSet_figures/", brain_region,
                                         enrich_results, DE_result, up_or_down, ...){
  # select genesets from an enrichment result by key words
  # extract member genes' ID from selected genesets
  # extract the FC of these genes from a DE result (for bar plot)
  # also return a list of DEGs related to key words
  
  dir.create(fig_dir)
  
  mem_DEGs_all <- c()
  gene_sets_all <- c()
  gs_table_all <- tibble(gene_set = character(), member_DEG = character())
  
  for(enrich_result in enrich_results){  # e.g. test_GSA$reactome@result
    mem_DEGs <- enrich_result@result[str_detect(enrich_result@result$Description,
                                                regex(key_word, ignore_case = TRUE)), 
                                     "geneID"]
    gene_sets <- enrich_result@result[str_detect(enrich_result@result$Description,
                                                 regex(key_word, ignore_case = TRUE)), 
                                     "Description"]
    gs_table <- tibble(gene_set = gene_sets, member_DEG = mem_DEGs)
    if(dim(gs_table)[1] != 0){
      for(i in 1:dim(gs_table)[1]){
        gst_tmp <- strsplit(gs_table$member_DEG[i], split = "/")[[1]]
        gst_tmp <- filter(symbol_entrezid, entrezid %in% gst_tmp)$symbol
        gs_table$member_DEG[i] <- paste0(gst_tmp, collapse = ", ")
      }
      gs_table_all <- bind_rows(gs_table_all, gs_table)  # cannot combine empty tibble to empty tibble
    }
    
      
    mem_DEGs <- strsplit(paste0(mem_DEGs, collapse = "/"), split = "/")[[1]]
    #mem_DEGs <- unique.default(mem_DEGs)
    mem_DEGs_all <- unique.default(c(mem_DEGs_all, mem_DEGs))
    gene_sets_all <- unique.default(c(gene_sets_all, gene_sets))
  }
  
  
  mem_DEGs_sym <- filter(symbol_entrezid, entrezid %in% mem_DEGs_all)$symbol
  
  gs_table_all$gene_set <- gsub("_", " ", gs_table_all$gene_set)
  # too long for table cells if keeping the "_"
  
  # save relevant DEGs and the geneset - DEG table locally
  saveRDS(mem_DEGs_sym, file = file.path(fig_dir, paste0(brain_region, "_", 
                                                         fig_name, "_DEGs.rds")))
  saveRDS(gs_table_all, file = file.path(fig_dir, paste0(brain_region, "_", 
                                                         fig_name, "_geneset_table.rds")))
  
  # also save the geneset - DEG table to LaTeX format
  print(xtable(gs_table_all, type = "latex", align = "cp{9.5cm}p{6.5cm}",
               caption = paste0("genesets related to ", fig_name, 
                                " and their member DEGs", 
                                ", ", gsub("_", "\\_", brain_region, fixed = TRUE))), 
                                # even "_" need to be protected by "\" in LaTeX
        file = file.path(fig_dir, paste0(brain_region, "_", fig_name, 
                                         "_geneset_table.tex")),
        booktabs = TRUE, tabular.environment = "longtable", width = "\\textwidth",
        caption.placement = "top")  
        # align = "cXX"; tabular.environment = tabularx, table.placement = "H"
  
  # prepare for plotting the error bar
  if(up_or_down == "up"){
    idx <- DE_result$log2FoldChange > 0
  }
  else if(up_or_down == "down"){
    idx <- DE_result$log2FoldChange < 0
  }
  idx[is.na(idx)] <- FALSE
  DE_result <- DE_result[idx, ]
  
  DE_result <- as_tibble(DE_result, rownames = "symbol")
  #DE_result$FC <- 2 ** DE_result$log2FoldChange
  DE_result <- filter(DE_result, symbol %in% mem_DEGs_sym)
  DE_result <- inner_join(symbol_entrezid, DE_result, by = "symbol")
  DE_result[DE_result[, "log2FoldChange"] < 0, "lfcSE"] <- DE_result[DE_result[, "log2FoldChange"] < 0, "lfcSE"] * (-1)
  
  ### plotting using ggplot2
  
  # bar plot using DE_result
  ggplot(data = DE_result) + 
    #geom_bar(mapping = aes(x = symbol, y = FC - 1), stat = "identity")
    geom_bar(mapping = aes(y = reorder(symbol, abs(log2FoldChange)), 
                           x = log2FoldChange), 
             stat = "identity") +
    # theme(axis.text.y = element_text(angle = 90, hjust = 0.3)) +
    geom_errorbar(mapping = aes(y = reorder(symbol, abs(log2FoldChange)), 
                                x = log2FoldChange,
                                xmin = log2FoldChange, xmax = log2FoldChange + lfcSE), 
                  width=0.8) + 
    ylab("Gene") +
    theme(axis.text.y = element_text(size = 7)) + 
    ggsave(path = fig_dir, units = "in", width = 7, ..., # pass the "..." parameter here
           paste0(brain_region, "_", fig_name, "_log2fc.png"))
  # use Log2FC to plot, and use LfcSE to add the error bar

}


CompactGeneTable <- function(fig_name, fig_dir = "GeneSet_figures/", brain_region){
  # as the geneset - member DEG tables are only suitable for supplementary table,
  # this function makes a up - down x male - female table for each key word
  # also in tex format
  # this time the brain region does not include gender info
  
  gene_table <- tibble(gender = c("male", "female"), 
                       up = c("", ""), down = c("", ""))
  
  M_up <- get(paste0(brain_region, "_M", "_up_", fig_name, "_geneset_table"))
  M_down <- get(paste0(brain_region, "_M", "_down_", fig_name, "_geneset_table"))
  F_up <- get(paste0(brain_region, "_F", "_up_", fig_name, "_geneset_table"))
  F_down <- get(paste0(brain_region, "_F", "_down_", fig_name, "_geneset_table"))
  
  
  gene_table$up[1] <- paste0(sort(unique.default(strsplit(paste0(M_up$member_DEG, 
                                                     collapse = ", "),
                                              split = ", ")[[1]])), collapse = ", ")
  gene_table$down[1] <- paste0(sort(unique.default(strsplit(paste0(M_down$member_DEG, 
                                                       collapse = ", "),
                                                split = ", ")[[1]])), collapse = ", ")
  gene_table$up[2] <- paste0(sort(unique.default(strsplit(paste0(F_up$member_DEG, 
                                                     collapse = ", "),
                                              split = ", ")[[1]])), collapse = ", ")
  gene_table$down[2] <- paste0(sort(unique.default(strsplit(paste0(F_down$member_DEG, 
                                                       collapse = ", "),
                                                split = ", ")[[1]])), collapse = ", ")
  
  saveRDS(gene_table, file = file.path(fig_dir, paste0(brain_region, "_", fig_name, 
                                                       "_compact_table.rds")))
  
  print(xtable(gene_table, type = "latex", align = "p{0cm}p{2cm}p{7cm}p{7cm}",
               caption = paste0("DEGs related to ", fig_name, 
                                ", ", gsub("_", "\\_", brain_region, fixed = TRUE))), 
        # even "_" need to be protected by "\" in LaTeX
        file = file.path(fig_dir, paste0(brain_region, "_", fig_name, 
                                         "_compact_table.tex")),
        booktabs = TRUE, tabular.environment = "longtable",
        caption.placement = "top")
  
}





VolcanoMAplotGenes <- function(key_genes, DE_result, 
                               fig_dir = "GeneSet_figures/", 
                               brain_region, fig_name, 
                               alpha = 0.05){
  # make MA plot and volcano plot highlighting important DEGs
  
  
  
  # MA plot
  png(file = paste0(fig_dir, "MA_", brain_region, "_", fig_name, ".png"),
      units = "in", width = 7, height = 7, res = 300)
  DESeq2::plotMA(DE_result, alpha = alpha, main = "", ylim=c(-1, 1))
  dev.off()
}



# apply the functions to each brain regions ---------------------------

setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/")

##### frontal lobe (FL)


### load the dds and res objects for male and female
dds_combined_FL_F_batches_sva_vst <- readRDS("DE_combined_batches_sva_figures/dds_combined_FL_F_batches_sva_vst.rds")
dds_combined_FL_M_batches_sva_vst <- readRDS("DE_combined_batches_sva_figures/dds_combined_FL_M_batches_sva_vst.rds")
DEres_combined_FL_F_braak_batches_sva_vst <- readRDS("DE_combined_batches_sva_figures/DEres_combined_FL_F_braak_batches_sva_vst.rds")
DEres_combined_FL_M_braak_batches_sva_vst <- readRDS("DE_combined_batches_sva_figures/DEres_combined_FL_M_braak_batches_sva_vst.rds")

### get the background genes (genes that are actually expressed in the certain brain regions of certain gender)

bg_gene_M <- DetermineBGgenes(dds_combined_FL_M_batches_sva_vst, symbol_entrezid)
bg_gene_F <- DetermineBGgenes(dds_combined_FL_F_batches_sva_vst, symbol_entrezid)


### geneset analysis 


# default alpha = 0.05, lfc_threshold = 0.1 (for absolute values)
GeneSetAnalysis(res = DEres_combined_FL_M_braak_batches_sva_vst,  # 3min
                brain_region = "combined_FL_M",
                bg_gene = bg_gene_M, wpid2gene, wpid2name, ms2gene)

GeneSetAnalysis(res = DEres_combined_FL_F_braak_batches_sva_vst,  
                brain_region = "combined_FL_F",
                bg_gene = bg_gene_F, wpid2gene, wpid2name, ms2gene)


### bar plot & geneset table for each keywords, each sex and regulation direction

geneset_dir = "GeneSet_figures"
# load the results of geneset analysis
for (geneset_list in list.files(path = geneset_dir, pattern = "geneset_list")) {
  assign(gsub("(^.+)\\.rds", "\\1", geneset_list), 
         readRDS(file.path(geneset_dir, geneset_list)))
}

keyWordsTable <- tribble(~key_word,                           ~fig_name,   ~heightUp, ~heightDown,
                         "mitochondri",                       "mitochon",  7,         20,
                         "mitophag",                          "mitophagy", 3,         7,
                         "lysoso",                            "lysosome",  7,         7,
                         "oxidative stress|oxidative_stress", "oxid",      8,         8)


for (gsl in ls(pattern = "combined.+geneset_list.+")){  # M and F are separated
  sex_info <- gsub("combined_.+_([F|M])_geneset_list_.+", "\\1", gsl)
  dataset_info <- gsub("(combined_.+)_[F|M]_geneset_list_.+", "\\1", gsl)
  up_or_down <- gsub("combined_.+_[F|M]_geneset_list_(.+)", "\\1", gsl)
  
  for (i in 1:dim(keyWordsTable)[1]) {
    VisualiseGeneChange_keywords(key_word = keyWordsTable$key_word[i], 
                                 fig_name = keyWordsTable$fig_name[i], 
                                 up_or_down = up_or_down,
                                 symbol_entrezid = symbol_entrezid, 
                                 enrich_results = get(paste0(dataset_info, "_",
                                                             sex_info, "_",
                                                             "geneset_list_",
                                                             up_or_down)),
                                 DE_result = get(paste0("DEres_",
                                                        dataset_info, "_",
                                                        sex_info, "_",
                                                        "braak_batches_sva_vst")),
                                 brain_region = paste0(dataset_info, "_", 
                                                       sex_info, "_",
                                                       up_or_down),
                                 height = if_else(up_or_down == "up", 
                                                  keyWordsTable$heightUp[i], 
                                                  keyWordsTable$heightDown[i]))  
                                 # pass to ggsave as "..."
  }
  
}


system("mkdir GeneSet_figures/genesets")
system("mv GeneSet_figures/*_log2fc.png GeneSet_figures/genesets")
system("mv GeneSet_figures/*_geneset_table.tex GeneSet_figures/genesets")



### compact gene table for each key word

geneset_dir = "GeneSet_figures"
# load the geneset - member DEG table
for (geneset_table in list.files(path = geneset_dir, pattern = "geneset_table")) {
  assign(gsub("(^.+)\\.rds", "\\1", geneset_table), 
         readRDS(file.path(geneset_dir, geneset_table)))
}


CompactGeneTable(fig_name = "lysosome", fig_dir = "GeneSet_figures/", 
                 brain_region = "combined_FL")

CompactGeneTable(fig_name = "oxid", fig_dir = "GeneSet_figures/", 
                 brain_region = "combined_FL")

CompactGeneTable(fig_name = "mitophagy", fig_dir = "GeneSet_figures/", 
                 brain_region = "combined_FL")

CompactGeneTable(fig_name = "mitochon", fig_dir = "GeneSet_figures/", 
                 brain_region = "combined_FL")


system("mv GeneSet_figures/*_compact_table.tex GeneSet_figures/genesets")


# for writing the Results section

combined_FL_mitochon_compact_table <- readRDS("GeneSet_figures/combined_FL_mitochon_compact_table.rds")
combined_FL_mitophagy_compact_table <- readRDS("GeneSet_figures/combined_FL_mitophagy_compact_table.rds")
combined_FL_oxid_compact_table <- readRDS("GeneSet_figures/combined_FL_oxid_compact_table.rds")
combined_FL_lysosome_compact_table <- readRDS("GeneSet_figures/combined_FL_lysosome_compact_table.rds")


CountAndAhowOverlap <- function(compact_table){
  # count and show overlapping genes in male and female; up- and down-regulated
  
  M_up <- strsplit(compact_table$up[1], split = ", ")[[1]]
  F_up <- strsplit(compact_table$up[2], split = ", ")[[1]]
  M_down <- strsplit(compact_table$down[1], split = ", ")[[1]]
  F_down <- strsplit(compact_table$down[2], split = ", ")[[1]]
  
  print(length(M_up))
  print(length(M_down))
  print(length(F_up))
  print(length(F_down))
  
  print(length(M_up[M_up %in% F_up]))
  print(paste0(M_up[M_up %in% F_up], collapse = ", "))
  print(length(M_down[M_down %in% F_down]))
  print(paste0(M_down[M_down %in% F_down], collapse = ", "))
}

CountAndAhowOverlap(combined_FL_lysosome_compact_table)




### MA plot and volcano plot

res = DEres_combined_FL_M_braak_batches_sva_vst

res <- filter(as_tibble(res, rownames = "gene"), padj < 0.05)
res <- filter(res, abs(log2FoldChange) >= 1)
key_genes <- res$gene

VolcanoMAplotGenes(key_genes = key_genes,  # key_genes above
                   DE_result = DEres_combined_FL_M_braak_batches_sva_vst,
                   fig_dir = "GeneSet_figures/", 
                   brain_region = "combined_FL_M", fig_name = "abs_log2FC_1", 
                   alpha = 0.05)



res = DEres_combined_FL_F_braak_batches_sva_vst

res <- filter(as_tibble(res, rownames = "gene"), padj < 0.05)
res <- filter(res, abs(log2FoldChange) >= 1)
key_genes <- res$gene

VolcanoMAplotGenes(key_genes = key_genes,  # key_genes above
                   DE_result = DEres_combined_FL_F_braak_batches_sva_vst,
                   fig_dir = "GeneSet_figures/", 
                   brain_region = "combined_FL_F", fig_name = "abs_log2FC_1", 
                   alpha = 0.05)




# have to draw volcano plot outside the function, as EnhancedVolcano() produce no plots if called inside a function


# male
res = DEres_combined_FL_M_braak_batches_sva_vst
res <- filter(as_tibble(res, rownames = "gene"), padj < 0.05)
res <- filter(res, abs(log2FoldChange) >= 1)
key_genes <- res$gene

DE_result = DEres_combined_FL_M_braak_batches_sva_vst
brain_region = "combined_FL_M"
fig_name = "abs_log2FC_1"
fig_dir = "GeneSet_figures/"
alpha = 0.05

gene_lab <- row.names(DE_result)
gene_lab[!(gene_lab %in% key_genes)] <- ""

png(file = paste0(fig_dir, "volano_", brain_region, "_", fig_name, ".png"),
    units = "in", width = 7, height = 7, res = 300)
EnhancedVolcano(DE_result, lab = gene_lab, x = 'log2FoldChange',
                y = 'padj', pCutoff = alpha, FCcutoff = 1.0,
                ylab = bquote(~-Log[10]~italic(padj)), 
                ylim = c(0, max(-log10(DE_result[['padj']]), na.rm=TRUE) + 0.1),
                xlim = c(min(DE_result[['log2FoldChange']], na.rm=TRUE) - 0.1,
                         max(DE_result[['log2FoldChange']], na.rm=TRUE) + 0.1),
                title = "", subtitle = "", pointSize = 0.1,
                legendLabels = c('NS', expression(Log[2]~FC),
                                 'padj', expression(padj~and~log[2]~FC)))
dev.off()


# female
res = DEres_combined_FL_F_braak_batches_sva_vst
res <- filter(as_tibble(res, rownames = "gene"), padj < 0.05)
res <- filter(res, abs(log2FoldChange) >= 1)
key_genes <- res$gene

DE_result = DEres_combined_FL_F_braak_batches_sva_vst
brain_region = "combined_FL_F"
fig_name = "abs_log2FC_1"
fig_dir = "GeneSet_figures/"
alpha = 0.05

gene_lab <- row.names(DE_result)
gene_lab[!(gene_lab %in% key_genes)] <- ""

png(file = paste0(fig_dir, "volano_", brain_region, "_", fig_name, ".png"),
    units = "in", width = 7, height = 7, res = 300)
EnhancedVolcano(DE_result, lab = gene_lab, x = 'log2FoldChange',
                y = 'padj', pCutoff = alpha, FCcutoff = 1.0,
                ylab = bquote(~-Log[10]~italic(padj)), 
                ylim = c(0, max(-log10(DE_result[['padj']]), na.rm=TRUE) + 0.1),
                xlim = c(min(DE_result[['log2FoldChange']], na.rm=TRUE) - 0.1,
                         max(DE_result[['log2FoldChange']], na.rm=TRUE) + 0.1),
                title = "", subtitle = "", pointSize = 0.1,
                legendLabels = c('NS', expression(Log[2]~FC),
                                 'padj', expression(padj~and~log[2]~FC)))
dev.off()

