# ScatterBoxplot_PubMed.R



# define function ---------------------------------------------------------



ScatterplotCertainGenes <- function(dds, gene_syms = c("ATG5", "OPTN", "PINK1"), 
                                    brain_region,
                                    fig_dir = "GeneSet_figures/"){
  
  mod <- model.matrix(design(dds), colData(dds))
  mat <- counts(dds, normalized = TRUE)
  mat <- limma::removeBatchEffect(x = mat,
                                  covariates = mod[, c(-1, -dim(mod)[2])])
  mat <- round(mat)
  mat[mat < 0] <- 0
  mode(mat) <- "integer"
  counts(dds) <- mat  # only for plotting, not suitable for DE analysis
  
  normalisedCounts <- tibble(count = numeric(), condition = character(), gene = character())
  for (gene in gene_syms) {
    normalisedCounts_tmp <- DESeq2::plotCounts(dds, gene = gene, intgroup = "condition",
                                               returnData = TRUE, normalized = FALSE,
                                               transform = FALSE)
    normalisedCounts_tmp$gene <- gene
    row.names(normalisedCounts_tmp) <- NULL
    normalisedCounts <- bind_rows(normalisedCounts, normalisedCounts_tmp)
  }
  
  normalisedCounts$condition <- factor(normalisedCounts$condition, levels = c("ctrl", "AD"))
  ggplot(data = normalisedCounts) +
    geom_boxplot(mapping = aes(x = condition, y = count)) +
    stat_compare_means(mapping = aes(x = condition, y = count),
                       method = "t.test") +
    # geom_signif(comparisons = list(c("ctrl", "AD")), 
    #             map_signif_level = TRUE) +  # does not work
    # default test = "wilcox.test", can also be t.test
    geom_point(mapping = aes(x = condition, y = count), 
               position = position_jitter(w = 0.1, h = 0), 
               size = 1) +
    facet_wrap(~ gene, scales = "free") + 
    ggsave(paste0(fig_dir, brain_region, "_counts_", 
                  paste0(gene_syms, collapse = "_"), ".png"),
           height = 7, width = 7, units = "in")
  
}




# apply functions ---------------------------------------------------------


setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna")

dir.create("Scatter_Boxplot")

### lysosome

## overlap between male and female PPI 

# VAMP2, 15 publications
# member of Reactome lysosome vesicle biogenesis pathway
# M: 0.00028
ScatterplotCertainGenes(dds = dds_combined_FL_M_batches_sva_vst, 
                        gene_syms = c("VAMP2"), 
                        brain_region = "combined_FL_M",
                        fig_dir = "Scatter_Boxplot/")

# F: 0.014
ScatterplotCertainGenes(dds = dds_combined_FL_F_batches_sva_vst, 
                        gene_syms = c("VAMP2"),
                        brain_region = "combined_FL_F",
                        fig_dir = "Scatter_Boxplot/")

## male PPI 

# M6PR, only 2 recent publications,
# "Importantly, we have also demonstrated increased alternative splicing of APMAP 
# and lowered levels of the Aβ controllers HSPA1A and CD-M6PR 
# in human brains from neuropathologically verified AD cases."

# member of Reactome lysosome vesicle biogenesis pathway
# M: 0.00016
ScatterplotCertainGenes(dds = dds_combined_FL_M_batches_sva_vst, 
                        gene_syms = c("M6PR"), 
                        brain_region = "combined_FL_M",
                        fig_dir = "Scatter_Boxplot/")

# F: N.S.
ScatterplotCertainGenes(dds = dds_combined_FL_F_batches_sva_vst, 
                        gene_syms = c("M6PR"),
                        brain_region = "combined_FL_F",
                        fig_dir = "Scatter_Boxplot/")




### mitophagy
# By STRING

## overlap between male and female

# GAPDH, member of all PPI and KEGG AD pathway, 89 papers
# M: N.S.
ScatterplotCertainGenes(dds = dds_combined_FL_M_batches_sva_vst, 
                        gene_syms = c("GAPDH"), 
                        brain_region = "combined_FL_M",
                        fig_dir = "Scatter_Boxplot/")

# F: 0.014
ScatterplotCertainGenes(dds = dds_combined_FL_F_batches_sva_vst, 
                        gene_syms = c("GAPDH"),
                        brain_region = "combined_FL_F",
                        fig_dir = "Scatter_Boxplot/")


# MITF, 4 papers
# M: 0.021
ScatterplotCertainGenes(dds = dds_combined_FL_M_batches_sva_vst, 
                        gene_syms = c("MITF"), 
                        brain_region = "combined_FL_M",
                        fig_dir = "Scatter_Boxplot/")

# F: 0.021
ScatterplotCertainGenes(dds = dds_combined_FL_F_batches_sva_vst, 
                        gene_syms = c("MITF"),
                        brain_region = "combined_FL_F",
                        fig_dir = "Scatter_Boxplot/")



### oxidative stress response

# PRDX2, 12 publication
# member of the Reactome deregulated CDK5 triggers multiple 
# neurodegenerative pathways in Alzheimer's disease models pathway
# M: 0.01
ScatterplotCertainGenes(dds = dds_combined_FL_M_batches_sva_vst, 
                        gene_syms = c("PRDX2"), 
                        brain_region = "combined_FL_M",
                        fig_dir = "Scatter_Boxplot/")

# F: 0.00056
ScatterplotCertainGenes(dds = dds_combined_FL_F_batches_sva_vst, 
                        gene_syms = c("PRDX2"),
                        brain_region = "combined_FL_F",
                        fig_dir = "Scatter_Boxplot/")


