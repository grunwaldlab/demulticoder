library(devtools)
#if you need to make changes to functions in R sub-directory, and then use load_all command again, and then re-document with document()
setwd("~/rps10_metabarcoding_tool/test_package") #be careful about this
load_all("~/rps10_metabarcoding_tool/test_package")
document()

#Example
#Note-sample names should be sample1_xx_xx_R1.fastq.gz. The key is all sample names are unique before first underscore, and this sample name matches metadata.csv sample name.

#You will need to modify directory names and file names
directory_path<-"~/WindRiver_test_mas/" ##choose a directory for all downstream steps
raw_path <-file.path(directory_path, "raw_data") #place your raw data files, csv files, and downloaded databases into raw_data subdirectory into your main directory
primer_path <-file.path(raw_path, "primers.csv") ##modify .csv name or keep this name
#Metadata file just needs sample_name one column, and primer_name in second column (this function is being tweaked-see example)
metadata_path <-file.path(raw_path, "metadata.csv") ##modify .csv name or keep this name. The sample_name in the metadata sheet needs to match the first part (before first underscore), of the zipped raw FASTQ files
cutadapt_path<-"/Users/masudermann/miniconda3/bin/cutadapt"
intermediate_path <- create_intermediate(directory_path)
#create intermediate data folder in working directory
#Now you are ready for actual commands

#Further simplify
returnList<-prepare(directory_path, primer_path, metadata_path, fastq_path, intermediate_path, maxN=0, multithread=TRUE)
cut_trim(returnList,intermediate_path, cutadapt_path, verbose=TRUE, maxEE=5) 
#add message to let user know which steps have already run if they are skipped, and also provide info on the params. Make it easier for user to re-run analysis. 
asv_abund_matrix <- make_asvAbund_matrix(returnList, intermediate_path, returnList$cutadapt_data)

#Command to prepare databases for downstream steps
#TODO format like SILVA? 
create_ref_database(intermediate_path)
format_database(raw_path, "oomycetedb.fasta")
#If you included the ITS barcode in your analysis
#Should also test pipeline just on ITS dataset (could also make your own fungal database)
format_database2(raw_path, "unite.fasta")
#The remaining functions will be incorporated shortly

asv_pid_table<-process_rps10_barcode(returnList, asv_abund_matrix)

#Rename anything that is ambiguous