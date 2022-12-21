library(devtools)
#if you need to make changes to functions in R sub-directory, and then use load_all command again, and then re-document with document()
setwd("~/rps10_metabarcoding_tool/test_package") #be careful about this
load_all("~/rps10_metabarcoding_tool/test_package")
document()


#Example
#Note-sample names should be sample1_xx_xx_R1.fastq.gz. The key is all sample names are unique before first underscore, and this sample name matches metadata.csv sample name.

#You will need to modify directory names and file names
directory_path<-"~/Wen2" ##choose a directory for all downstream steps
raw_path <-file.path(directory_path, "raw_data") #place your raw data files, csv files, and downloaded databases into your chosen directory
primer_path <-file.path(raw_path, "primer_info.csv") ##modify .csv name or keep this name
#Metadata file just needs sample_name one column, and primer_name in second column (this function is being tweaked-see example)
metadata_path <-file.path(raw_path, "metadata5.csv") ##modify .csv name or keep this name. The sample_name in the metadata sheet needs to match the first part (before first underscore), of the zipped raw FASTQ files
cutadapt_path<-"/Users/masudermann/miniconda3/bin/cutadapt"
intermediate_path <- create_intermediate(directory_path)
#create intermediate data folder in working directory
#Now you are ready for actual commands


#You will need to modify directory names and file names
directory_path<-"~/rps10_metabarcoding_tool/test_package/data" ##choose a directory for all downstream steps
raw_path <-file.path(directory_path) #place your raw data files, csv files, and downloaded databases into your chosen directory
primer_path <-file.path(raw_path, "primer_info.csv") ##modify .csv name or keep this name
#Metadata file just needs sample_name one column, and primer_name in second column (this function is being tweaked-see example)
metadata_path <-file.path(raw_path, "metadata.csv") ##modify .csv name or keep this name. The sample_name in the metadata sheet needs to match the first part (before first underscore), of the zipped raw FASTQ files
cutadapt_path<-"/Users/masudermann/miniconda3/bin/cutadapt"
intermediate_path <- create_intermediate(directory_path)
#create intermediate data folder in working directory
#Now you are ready for actual commands
 
#Further simplify
#change names
returnList<-prepare_reads(directory_path, primer_path, metadata_path, fastq_path, intermediate_path, maxN=0, multithread=TRUE)
cut_trim(returnList,intermediate_path, cutadapt_path, verbose=TRUE, maxEE=5, truncQ=5, minLen=200, min_length=50) #fix cutadapt param!!!
#add message to let user know which steps have already run if they are skipped, and also provide info on the params. Make it easier for user to re-run analysis. 
#minboot? 
asv_abund_matrix <- make_asv_abund_matrix(returnList, intermediate_path, returnList$cutadapt_data, minOverlap=15, maxMismatch=2, verbose=TRUE, multithread=TRUE) #check when 0 again
#TODO adjust maxmismatch if there aren't alot of merged reads                                                                                         

#Wrapper function to format databases and assign taxonomy
summary<-assignTax_function(raw_path, "oomycetedb.fasta", "unite.fasta", returnList, asv_abund_matrix, multithread = TRUE, is_ITS=TRUE)
#output is weird
sessioninfo::session_info()

#outputs
