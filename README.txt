1. Quality control & expression quantification

QualityControlForGEO.R
    • apply fastqc on each datasets
	
test_trim_quant.sh
    • apply trimmomatic and salmon quant on each PE or SE datasets

2. DE analysis & preparation

Tximport.R
    • a function txiForEachBrainRegion() to apply tximport to each GEO and AMP-AD dataset, using transcripts – gene symbol annotation table from GENCODE v33;
    • apply the function to each dataset, save the txi object of each dataset as rds in rna_seq_for_mrna\;

colData_metadata.R
    • use hand-written GSE – SRP, SRP – SRR table from SRA_Run_Members.tab to make a list of SRR to be downloaded which are saved 
	  in to_be_downloaded_from_GEO.rds and used by DownloadConvertGEO.R, ContinueDownloadConvertGEO.R & ContinueDownloadConvertGEO_v2.R);
    • use GEOquery::getGEO() to get the ExpressionSet (eset) object for each GSE dataset, and save in esets_from_geo.rds;
    • examine the phenotype (gender, age, braak) of patients in each GSE dataset as is stored in pData and count patients to be included; 
	  do the same to AMP-AD datasets using metadata tables from Synapse;
    • use to_be_downloaded_from_GEO.txt to make a SRR – SRX – SAMN table;
    • a function add_SAMN_SRR_GSE() to add SAMN, SRR, GSE to each line of colData, and save the colData for the GSE as .rds file;
    • for each GSE, extract preliminary colData from pData of the eset, filter samples by age and braak, and use the function add_SAMN_SRR_GSE() to finalise it;
    • for each AMP-AD dataset, use metadata tables from Synapse to make colData by filtering samples by age and braak; save the colData as .rds file;

predict_gender.R
    • make transcript – gene symbol table from GENCODE v33	 annotation file gencode.v33.metadata.HGNC;
    • a function ChromAnno() to extract genes on chromosome X and/or Y from gencode.v33.annotation using the above mentioned table; apply the function;
    • a function OverviewXY() to make 2 DESeqDataSet (dds) objects using txi, colData and genes on chromosome X & Y respectively; 
	  make a PCA plot of the expression of genes on chromosome Y to see if the known gender annotation of samples correspond to 2 clusters 
	  (plots are saved in the folder rna_seq_for_mrna/); 
	  for each sample, calculate the total read counts of genes on chromosome Y divided by total read counts of genes on chromosome X;
    • a function DecideGender()  using above produced Y - X ratio and PCA data to predict the gender of each sample; samples in the cluster with a 
	  higher Y - X ratio are re-annotated as male, which is shown in a new column in colData; 
    • for each dataset, load existing txi and colData, and apply the 2 functions above to it; save the new colData containing inferred sex column 
	  as .rds file in the folder rna_seq_for_mrna/ (old ones are moved to the sub folder old_files/);

BatchEffectCorrection.R
    • a function NumeriseCovariates() to convert the covariates (g.e. age & pmi) in colData to numeric mode (default is factor), 
	  update the design formula (initially ~ condition) with the given one, and remove samples with NAs in these covariates from the given dds;
    • a function SVAtoDDS() to calculate given number (1, 2, 3 or 4) of surrogate variables (SVs) for the given dds, and add these SVs to 
	  the design formula (initially ~ condition);
    • a function BatchesSVAtoDDS()  to calculate given number (1, 2, 3 or 4) of surrogate variables (SVs) for the given dds, and add these SVs 
	  to the design formula while keeping existing covariates (e.g. age & pmi)
    • a function NegativeControlGenesFromRes() to determine negative control genes for male and female respectively using given dds of male and female 
	  and DESeq results object (res) for male and female, preparing for RUVseq;
    • a function GetWsByRUVg() to calculate the values of 2 RUVseq’s W variables for each samples, using dds, res and negative control genes of male and female;
    • (choose between “FALSE” i.e. do not correct, “sva”, “batches” i.e. known covariates and batches, “batches_sva” and “RUVseq”)

DESeq2_functions.R
    • source the script “BatchEffectCorrection.R”;
    • a function ConstructDDSandVisualise() to construct dds object from txi and colData after correcting the condition (diagnosis) info 
	  using braak score, divide the dds object into male and female, correct for batch effect (BE), variance stabilisation transformation, 
	  remove batch effect from transformed dds, visualise the transformed dds using PCA, dendrogram & heatmap, and save an untransformed dds 
	  (with correction variables in the design formula) for male and female respectively as .rds files in the sub folder for corresponding types of correction 
	  (i.e. DE_combined_batches_figures, DE_combined_batches_sva_figures, DE_combined_sva_figures, DE_combined_test_figures);
    • a function DEanalysis() to perform DE analysis on the male dds and female dds, extract the result objects and save as .rds files in the sub folder 
	  for corresponding types of correction, and make MA plot and Volcano plot using the result with a given alpha value (default 0.05), 
	  which are also saved in the same folders;

DESeq2_EachDataSets.R
    • source the script “DESeq2_functions.R” and load no extra R packages
    • for each dataset, load the saved txi and colData objects, use the function ConstructDDSandVisualise() and DEanalysis() to construct dds objects 
	  and perform DE analysis with or without different types of batch effect coorection; 
    • the dds and res objects, heatmap, PCA plots, hierarchical clustering tree, MA plot and Volcano plots are saved in folders corresponding to the types of BE
      correction methods (one of DE_test_figures/, DE_sva_figures/, DE_batches_figures/, DE_batches_sva_figures/ & DE_ruvseq_figures/)
    • all dataset are processed with no BE correction (products saved in DE_test_figures/) and with only sva correction (products saved in DE_sva_figures/); 
    • GSE110731, ROSMAP (FL) and MSBB (FL)  are also processed with batches + sva (products saved in DE_batches_sva_figures/); 
	  batches in MSBB (FL) are combined into just 3 according to clusters in the PCA plot;
    • batches and RUVseq (use batches + sva to get negative control genes) codes have been written for GSE110731, ROSMAP (FL) and MSBB (FL), 
	  including saving negative control genes as .rds files, but only RUVseq for GSE110731 and batches for GSE110731 and MSBB (FL) have been actually run;
    • a function BarplotVennplot() to visualise and compare the number of DEGs in male and female of several datasets from the same brain region 
	  using bar plot and Venn plot; the res objects are loaded from a given folder and assigned to correspondingly named variables using assign(), 
	  and accessed by a for loop using ls(pattern) on cached variables (objects in the workspace);
    • for Frontal Lobe (FL) apply the function BarplotVennplot() on GSE110731, ROSMAP (FL) and MSBB (FL) from folder DE_batches_sva_figures/, 
	  and save the figures in Bar_Venn/;
    • for Frontal Lobe (FL) archived the figures in Bar_Venn/, and re-run BarplotVennplot() on GSE110731, ROSMAP (FL), MSBB (FL) and combined (FL) 
	  after loading res for male and female of combined (FL) from folder DE_combined_batches_sva_figures/;

DESeq2_EachBrainRegions.R
    • a function Combine_ColDatas_Txis() to combine a list of colData into one colData and a list of txi to one txi; also specify variables to be kept 
	  as columns of colData; save the colData object and txi object as .rds files;
    • load the txi and colData of 3 frontal lobe datesets: GSE110731, ROSMAP (FL) and MSBB (FL), combine MSBB’s batches into 3 batches according to PCA plot, 
	  set ROSMAP (FL) and GSE110731 as another 2 batches; unify the Post-mortem interval as in minutes;
    • call the function Combine_ColDatas_Txis() on the 3 datasets, creating a combined txi object and a combined colData object;
    • call functions ConstructDDSandVisualise() and Deanalysis() to construct dds and perform DE analysis; results & figures without BE correction, 
	  with batches + sva correction, only batches correction, and with only sva correction are saved in DE_combined_test_figures/, 
	  DE_combined_batches_sva_figures/, DE_combined_batches_figures/, DE_combined_sva_figures/ respectively;

3. downstream analysis on DEGs

GenesetAnalysis.R
    • download and load .gmt file for wikiPathways database, load MsigDb from R package, and create necessary tables for overrepresented analysis;
    • make symbol – entrezid table from GENCODE v33 annotation files;
    • a function DetermineBGgenes() to select all genes with an average count above 1 across samples as background genes
    • a function GeneSetAnalysis() to extract DEGs and log2FC that meet a certain threshold (default is padj < 0.05) from given res object, 
	  split the DEG list into up- and down- regulated gene lists, use 6 databases for overrepresented analysis (Hypergeometric distribution) 
	  and save the 2 resulting lists of enriched genesets as 2 .rds files in folder GeneSet_figures/;the threshold parameters for geneset are 
	  pvalueCutoff = 1, qvalueCutoff = 1 in order to identify all genesets that include at least one member DEG; 
    • a function VisualiseGeneChange_keywords() to extract all DEGs related to provided key words from a list of up or down regulated genesets 
	  (including genesets from several databases), extract log2FC and lfcSE from res object, and visualise as bar plot; the selected DEG list and 
	  geneset – member DEG table are also saved as .rds files in the folder GeneSet_figures/, so are and other 4 .rds that keep up- and down- regulated gene 
	  separated; tables are also exported as LaTeX code;
    • a function ScatterplotCertainGenes() to extract the normalised counts after BE correction for given genes from a given dds object and visualise 
	  as scatter plot and box plot in one figurefigures are saved in the folder GeneSet_figures/;
    • a function HeatmapGenes() to make a heatmap by subsetting the given dds by specified key genes; the figure is saved at the folder GeneSet_figures/, 
	  with proper setting to ensure that the png() figure quality is the same to that of pdf() or ggplot2;
    • a function VolcanoMAplotGenes() to make an MA plot and a volcano plot from a res object while highlighting the given keep genes in the plot; 
	  figure are saved in GeneSet_figures/ with similar parameter values for png() as above;
    • load the dds and res object for male and female parts of the combined frontal lobe dataset, get the background genes 
	  and perform overrepresentation geneset analysis for male and female respectively; use the function DetermineBGgenes() and GeneSetAnalysis() each twice;
    • automatically load the 4 saved geneset lists (for up- and down- regulated genesets of male and female frontal lobe) using regular expression; 
	  call the function VisualiseGeneChange_keywords() on these geneset lists using one of the 4 key words 
	  (mitochondria, mitophagy, lysosome and oxidative stress response) each time on a for loop;
    • informally try the ScatterplotCertainGenes() function on the male and female’s dds object with 1~3 selected genes (related to key words or with largest FC);
    • automatically load the 16 DEG lists (for for up- and down- regulated genesets of male and female frontal lobe related to the 4 key words) 
	  using regular expression; call the function HeatmapGenes() on them with 2 dds objects to produce 16 heatmaps which are saved in the folder GeneSet_figures/;

PathwayTopoAnalysisVisual.R
    • a function PathwayVisualise() to map DEGs in a provided res object to 4 KEGG pathways and move the plots to the folder Pathway_figures/;
    • load the dds and res objects of the combined frontal lobe dataset (with male and female separated) and call the function PathwayVisualise() on them;
    • load all pathway databases in the R package graphite, extract all pathways related to the 4 key words plus base excision repair, 
	  and make a tsv containing log2FC for cytoscape visualisation; 

NetworkAnalysis.R
    • use the STRING API to download the tsv file of Protein-Protein Interaction (PPI) network for a given list of key genes related to a key words; 
	  also include other genes in the res object that interact with at least one third of the key genes (through combining several hundreds other DEGs 
	  with the key genes each time, and use the STRING API to extract PPI table from which only outside – inside interaction pairs are kept);
	  
ScatterBoxplot_PubMed.R
    • visualise the normalised read counts of selected gene in boxplots and scatter plots;
	
cytoscape_commands.txt
    • load the STRING data in tsv format into Cytoscape
    • apply proper style for the network
	