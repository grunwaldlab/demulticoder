library("devtools")
library("roxygen2")
library("testthat")
library("knitr")
library("dada2")
#usethis::create_package("~/rAMP")

#if you need to make changes to functions in R sub-directory, and then use load_all command again, and then re-document with document()
setwd("~/rps10_metabarcoding_tool") #be careful about this
load_all("~/rps10_metabarcoding_tool")
document()

#Simplify inputs
directory_path<-"~/rps10_metabarcoding_tool/data/rps10_ITS/" ##choose a directory for all downstream steps
primer_path <-file.path(directory_path, "primer_info.csv") ##modify .csv name or keep this name
#Metadata file just needs sample_name one column, and primer_name in second column (this function is being tweaked-see example)
metadata_path <-file.path(directory_path,"metadata.csv") ##modify .csv name or keep this name. The sample_name in the metadata sheet needs to match the first part (before first underscore), of the zipped raw FASTQ files
cutadapt_path<-"/Users/masudermann/miniconda3/bin/cutadapt"

#Further simplify
#change names


data_tables<-prepare_reads(directory_path, primer_path, metadata_path, maxN=0, multithread=TRUE)
cut_trim(data_tables,directory_path, cutadapt_path, verbose=TRUE, maxEE=8, truncQ=2, minLen=200, maxLen=297, minCutadaptlength=50) #fix cutadapt param!!!
#add message to let user know which steps have already run if they are skipped, and also provide info on the params. Make it easier for user to re-run analysis.
#minboot?
asv_abund_matrix <- make_asv_abund_matrix(data_tables, directory_path, data_tables$cutadapt_data, minOverlap=15, maxMismatch=2, verbose=TRUE, multithread=TRUE) #check when 0 again
#TODO adjust maxmismatch if there aren't alot of merged reads

#Wrapper function to format databases and assign taxonomy
#Change barcode to 'rps10', 'its', or 'rps10_its'
summary<-assignTax(directory_path, data_tables, asv_abund_matrix, multithread = TRUE, barcode = "rps10_its")

#To do-Martha
#Put example commands into user-friendly Rmarkdown document
#Generalize for work with with 16S barcode-format 16S database
#Try out test datasets
#compile outputs as example

#Hung
#Check if primers are still present on ASV and remove sequences that do
#Simplify and further generalize functions
#Get package on cluster
sessioninfo::session_info()
