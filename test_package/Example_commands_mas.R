library(devtools)
#if you need to make changes to functions in R sub-directory, and then use load_all command again, and then re-document with document()
load_all("~/rps10_metabarcoding_tool/test_package")
document()

#Example
#In a directory of your choosing, make a subdirectory called "raw_data". My directory was called Wen_test
#Place your Paired-end-fastq reads, your metadata.csv and your primer_info.csv files into this folder. You can also place the raw databases you downloaded.
#Note-sample names should be sample1_xx_xx_R1.fastq.gz. The key is all sample names are unique before first underscore, and this sample name matches metadata.csv sample name.
#I may simplify this file structure
#Error handling and parameter optimization in process. To tweak some of the DADA2 parameters, you will need to make changes to the original functions in the R folder.

#You will need to modify directory names and file names
directory_path<-"~/WindRiver_test_mas" ##choose a directory for all downstream steps
raw_path <-file.path(directory_path, "raw_data") #place your raw data files, csv files, and downloaded databases to raw_data subdirectory into your main directory
primer_path <-file.path(raw_path, "primers_format.csv") ##modify .csv name or keep this name
#Metadata file just needs sample_name in first column, and primer_name in second column (this function is being tweaked-see example)
metadata_path <-file.path(raw_path, "metadata.csv") ##modify .csv name or keep this name
cutadapt_path<-"/Users/masudermann/miniconda3/bin/cutadapt"
intermediate_path <- create_intermediate(directory_path)
#create intermediate data folder in working directory
#Now you are ready for actual commands

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
ITS_rps10_barcode_function(intermediate_path, cutadapt_data)
#still some glitches with final steps

#TODO
#Parameters optimization
#Further clarify the inputs
#Templates for input files
#Consistency in naming
#Improve function names, and continue with documentation
#Check filter parameters-a lot of rps10 reads are being thrown out at the filter steps if you use default parameters
