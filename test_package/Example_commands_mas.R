library(devtools)
#if you need to make changes to functions in R sub-directory, and then use load_all command again, and then re-document with document()
setwd("~/rps10_metabarcoding_tool/test_package") #be careful about this
load_all("~/rps10_metabarcoding_tool/test_package")
document()

#Example
#Note-sample names should be sample1_xx_xx_R1.fastq.gz. The key is all sample names are unique before first underscore, and this sample name matches metadata.csv sample name.

#make interactive
directory_path<-"~/WindRiver_test_mas" ##choose a directory for all downstream steps
raw_path <-file.path(directory_path, "raw_data") #place your raw data files, csv files, and downloaded databases to raw_data subdirectory into your main directory
primer_path <-file.path(raw_path, "primers.csv") ##modify .csv name or keep this name
#Metadata file just needs sample_name in first column, and primer_name in second column (this function is being tweaked-see example)
metadata_path <-file.path(raw_path, "metadata.csv") ##modify .csv name or keep this name
#In a directory of your choosing, make a subdirectory called "raw_data". My directory was called Wen_test
#Place your Paired-end-fastq reads, your metadata.csv and your primer_info.csv files into this folder. You can also place the raw databases you downloaded.
#Note-sample names should be sample1_xx_xx_R1.fastq.gz. The key is all sample names are unique before first underscore, and this sample name matches metadata.csv sample name, which would just be in the spreadsheet as sample1, etc.
#I may simplify this file structure
#Error handling and parameter optimization in process. To tweak some of the DADA2 parameters, you will need to make changes to the original functions in the R folder.

#You will need to modify directory names and file names
directory_path<-"~/WindRiver_test_mas" ##choose a directory for all downstream steps
raw_path <-file.path(directory_path, "raw_data") #place your raw data files, csv files, and downloaded databases into raw_data subdirectory into your main directory
primer_path <-file.path(raw_path, "primers_format.csv") ##modify .csv name or keep this name
#Metadata file just needs sample_name one column, and primer_name in second column (this function is being tweaked-see example)
metadata_path <-file.path(raw_path, "metadata.csv") ##modify .csv name or keep this name. The sample_name in the metadata sheet needs to match the first part (before first underscore), of the zipped raw FASTQ files
cutadapt_path<-"/Users/masudermann/miniconda3/bin/cutadapt"
intermediate_path <- create_intermediate(directory_path)
#create intermediate data folder in working directory
#Now you are ready for actual commands

#Further simplify
returnList<-prepare(directory_path, primer_path, metadata_path, fastq_path, intermediate_path, maxN=0, multithread=TRUE)
cut_trim(returnList,intermediate_path, cutadapt_path, verbose=TRUE, maxEE=2) 
#add message to let user know which steps have already run if they are skipped, and also provide info on the params. Make it easier for user to re-run analysis. 

#Command to prepare databases for downstream steps
#TODO format like SILVA? 
create_ref_database(intermediate_path)
format_database(raw_path, "oomycetedb.fasta")
#If you included the ITS barcode in your analysis
#Should also test pipeline just on ITS dataset (could also make your own fungal database)
format_database2(raw_path, "unite.fasta")
#The remaining functions will be incorporated shortly.

asv_abund_matrix<-make_asvAbund_matrix(returnList, intermediate_path)

process_rps10_ITS_barcode <- function(returnList, intermediate_path, asv_abund_matrix)
#To run functions individually, run each of these functions below. You shouldn't have to change anything here. Note, for the DADA2 functions, we set some default parameters initially, but we are now working on generalizing this . Stay tuned. 
intermediate_path <- create_intermediate(directory_path)
primer_data <- prepare_primers(primer_path)
metadata <- prepare_metadata(metadata_path, primer_data)
fastq_data <- prepare_fastq(raw_path, intermediate_path)
pre_primer_hit_data<- pre_primer_hit_data(primer_data, fastq_data, intermediate_path)
pre_primer_plot <- primer_hit_plot(pre_primer_hit_data, fastq_data, intermediate_path, "pre_primer_plot.pdf")
cutadapt_data <- cutadapt_tibble(fastq_data, metadata, intermediate_path)
cutadapt_run(cutadapt_path, cutadapt_data)
quality_plots<-plot_qc(cutadapt_data, intermediate_path)
filter_results <-filter_and_trim(intermediate_path, cutadapt_data)
post_primer_hit_data <- post_primer_hit_data(primer_data, cutadapt_data, intermediate_path)
post_primer_plot <- primer_hit_plot(post_primer_hit_data, fastq_data, intermediate_path, "post_primer_plot.pdf")
quality_plots2 <- post_trim_qc(cutadapt_data, intermediate_path)
return(cutadapt_data)

#To run the first functions all together use this main function, and skip the individual functions above. 
cutadapt_data<-main_cutadapt_function(directory_path, primer_path, metadata_path, fastq_path,intermediate_path, cutadapt_path)

#Command to prepare databases for downstream steps
create_ref_database(intermediate_path)
format_database(raw_path, "oomycetedb.fasta")
#If you included the ITS barcode in your analysis
#Should also test pipeline just on ITS dataset
format_database2(raw_path, "unite.fasta")
#The remaining functions will be incorporated shortly.

#Command for just rps10 barcode, multiple samples
rps10_barcode_function(intermediate_path, cutadapt_data)
#TODO
#Command for rps10 and ITS barcodes, multiple samples
#This is still being revised. We are having some difficulities with running separate_abund_table function with Ricardo's dataset. It worked with another test dataset.
ITS_rps10_barcode_function(intermediate_path, cutadapt_data)

#TODO
#Parameters optimization
#Further clarify the inputs
#Templates for input files
#Consistency in naming
#Improve function names, and continue with documentation
#Check filter parameters-a lot of rps10 reads are being thrown out at the filter steps if you use default parameters
