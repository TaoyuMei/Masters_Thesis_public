# QualityConotrolForGEO.R
library(fastqcr)
library(ngsReports)
library(stringr)



# dealing with GSE raw data ------------------------------------------

## generating reports for each fastq files
setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna")
GSEs<- grep(pattern = "GSE.*", x = list.files(), value = TRUE)
for(gse in GSEs){
  fastqc(fastqc.path = "/binf-isilon/alab/students/vrw936/software/FastQC/fastqc",
         fq.dir = gse)
}

for (gse in c("GSE104704", "GSE110731", "GSE125050",
              "GSE125583", "GSE53697", "GSE95587")) {
  qc_report(qc.path = paste0("./", gse, "/FASTQC"),
            result.file = paste0("./", gse, "/FASTQC/", gse, "multi-qc-report"),
            interpret = TRUE, experiment = gse)

}



# fastqcr and ngsReports for the 6 GSE after trimmomatic ------------------

setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna")
GSEs<- grep(pattern = "GSE.*", x = list.files(), value = TRUE)

for(gse in GSEs){
  print(gse)
  fastqc(fastqc.path = "/binf-isilon/alab/students/vrw936/software/FastQC/fastqc",
         fq.dir = file.path(gse, "trimmed"))
}

altTemplate <- file.path("/binf-isilon", "alab", "students", "vrw936",
                         "Master_Thesis", "MyCode",
                         "ngsReports_Fastqc_template_no_overrepresent.Rmd")
for (gse in GSEs) {
  fileDir <- file.path(".", gse, "trimmed", "FASTQC")
  writeHtmlReport(fileDir, overwrite = TRUE, template = altTemplate)
}


# QC for each AMP-AD datasets ---------------------------------------------

### functionalise
qcAndReport <- function(wd = "/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna", 
                        fastq_dir){
  setwd(wd)
  fastqc(fastqc.path = "/binf-isilon/alab/students/vrw936/software/FastQC/fastqc", 
         fq.dir = fastq_dir)
}


### each dataset, before and after trimmomatic

# ROSMAP before trimmomatic
# moving ROSMAP fastqc results into 7 folders
for(i in 1:7){
 writeHtmlReport(file.path("ROSMAP/fastq/", "FASTQC", paste0("FASTQC_", i)),
                 overwrite = TRUE, template = altTemplate)
}


altTemplate <- file.path("/binf-isilon", "alab", "students", "vrw936",
                        "Master_Thesis", "MyCode",
                        "ngsReports_Fastqc_template_no_overrepresent.Rmd")

qcAndReport(fastq_dir = "MSBB/fastq/")
qcAndReport(fastq_dir = "MayoRNAseq/fastq/Mayo_CBE_sample_FASTQs/")
qcAndReport(fastq_dir = "MayoRNAseq/fastq/Mayo_TCX_sample_FASTQs/")



# ROSMAP after trimmomatic
qcAndReport(fastq_dir = "ROSMAP/fastq/trimmed")
setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/ROSMAP/fastq/trimmed/FASTQC")

for(i in 1:6){
 system(paste0("mkdir ", "FASTQC_", i))
 system(paste0("mv `ls | grep \".[html|zip]\" | head -360` ", "FASTQC_", i))
}
system("mkdir FASTQC_7")
system("mv `ls | grep \".[html|zip]\"` FASTQC_7")

altTemplate <- file.path("/binf-isilon", "alab", "students", "vrw936",
                        "Master_Thesis", "MyCode",
                        "ngsReports_Fastqc_template_no_overrepresent.Rmd")

for(i in 1:7){
 writeHtmlReport(file.path("ROSMAP/fastq/", "trimmed", "FASTQC", paste0("FASTQC_", i)),
                 overwrite = TRUE, template = altTemplate)
}


# functionalise
SeparateAndngsReports <- function(wd, num_fol){
  # a function to separate the FASTQC results into several folders
  # and combine the results using ngsReports
  # only for web RStudio in the server, not for background processes
  
  setwd(wd)
  for(i in 1:(num_fol - 1)){
    system(paste0("mkdir ", "FASTQC_", i))
    system(paste0("mv `ls | grep \".[html|zip]\" | head -360` ", "FASTQC_", i))
  }
  # actually 180 htmls and 180 zips for each folder, corresponding to 180 fastqs
  
  system(paste0("mkdir FASTQC_", num_fol))
  system(paste0("mv `ls | grep \".[html|zip]\"` FASTQC_", num_fol))
  
  altTemplate <- file.path("/binf-isilon", "alab", "students", "vrw936",
                           "Master_Thesis", "MyCode",
                           "ngsReports_Fastqc_template_no_overrepresent.Rmd")
  
  for(i in 1:num_fol){
    writeHtmlReport(file.path(paste0("FASTQC_", i)),
                    overwrite = TRUE, template = altTemplate)
  }
  
  
}


# MSBB before trimmmatic
SeparateAndngsReports(wd = "/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/MSBB/fastq/FASTQC",
                     num_fol = 6)

setwd("/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/MSBB/fastq/FASTQC")
for(i in 1:5){
 system(paste0("mkdir ", "FASTQC_", i))
 system(paste0("mv `ls | grep \".[html|zip]\" | head -360` ", "FASTQC_", i))
}

system("mkdir FASTQC_6")
system("mv `ls | grep \".[html|zip]\"` FASTQC_6")

altTemplate <- file.path("/binf-isilon", "alab", "students", "vrw936",
                        "Master_Thesis", "MyCode",
                        "ngsReports_Fastqc_template_no_overrepresent.Rmd")

for(i in 1:6){
 writeHtmlReport(file.path("MSBB/fastq/", "FASTQC", paste0("FASTQC_", i)),
                 overwrite = TRUE, template = altTemplate)
}



# MSBB after trimmomatic
qcAndReport(fastq_dir = "MSBB/fastq/trimmed")
SeparateAndngsReports(wd = "/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/MSBB/fastq/trimmed/FASTQC",
                     num_fol = 6)



# MayoRNAseq TCX before trimmomatic
SeparateAndngsReports(wd = "/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_TCX_sample_FASTQs/FASTQC",
                     num_fol = 3)


# MayoRNAseq TCX after trimmomatic
qcAndReport(fastq_dir = "MayoRNAseq/fastq/Mayo_TCX_sample_FASTQs/trimmed")
SeparateAndngsReports(wd = "/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_TCX_sample_FASTQs/trimmed/FASTQC",
                     num_fol = 3)


# MayoRNAseq CBE before trimmomatic
SeparateAndngsReports(wd = "/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_CBE_sample_FASTQs/FASTQC",
                     num_fol = 3)


# MayoRNAseq CBE after trimmomatic
qcAndReport(fastq_dir = "MayoRNAseq/fastq/Mayo_CBE_sample_FASTQs/trimmed")
SeparateAndngsReports(wd = "/binf-isilon/alab/students/vrw936/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_CBE_sample_FASTQs/trimmed/FASTQC",
                     num_fol = 3)



