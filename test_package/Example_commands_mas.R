library(devtools)
#if you need to make changes to functions in R sub-directory, and then use load_all command again, and then re-document with document()
setwd("~/rps10_metabarcoding_tool/test_package") #be careful about this
load_all("~/rps10_metabarcoding_tool/test_package")
document()

#Example
#In a directory of your choosing, make a subdirectory called "raw_data". My directory was called Wen_test
#Place your Paired-end-fastq reads, your metadata.csv and your primer_info.csv files into this folder. You can also place the raw databases you downloaded.
#Note-sample names should be sample1_xx_xx_R1.fastq.gz. The key is all sample names are unique before first underscore, and this sample name matches metadata.csv sample name.
#I may simplify this file structure
#Error handling and parameter optimization in process. To tweak some of the DADA2 parameters, you will need to make changes to the original functions in the R folder.

#You will need to modify directory names and file names
directory_path<-"~/Wen_test" ##choose a directory for all downstream steps
raw_path <-file.path(directory_path, "raw_data") #place your raw data files, csv files, and downloaded databases to raw_data subdirectory into your main directory
primer_path <-file.path(raw_path, "primer_info.csv") ##modify .csv name or keep this name
#Metadata file just needs sample_name in first column, and primer_name in second column (this function is being tweaked-see example)
metadata_path <-file.path(raw_path, "metadata5.csv") ##modify .csv name or keep this name
cutadapt_path<-"/Users/masudermann/miniconda3/bin/cutadapt"
intermediate_path <- create_intermediate(directory_path)
#create intermediate data folder in working directory
#Now you are ready for actual commands

directory_path<-"~/WindRiver_test_mas" ##choose a directory for all downstream steps
raw_path <-file.path(directory_path, "raw_data") #place your raw data files, csv files, and downloaded databases to raw_data subdirectory into your main directory
primer_path <-file.path(raw_path, "primers.csv") ##modify .csv name or keep this name
#Metadata file just needs sample_name in first column, and primer_name in second column (this function is being tweaked-see example)
metadata_path <-file.path(raw_path, "metadata.csv") ##modify .csv name or keep this name
cutadapt_path<-"/Users/masudermann/miniconda3/bin/cutadapt"
intermediate_path <- create_intermediate(directory_path)
#create intermediate data folder in working directory
#Now you are ready for actual commands


#each command individually
intermediate_path <- create_intermediate(directory_path)
primer_data <- prepare_primers(primer_path)
metadata <- prepare_metadata(metadata_path, primer_data)
fastq_data <- prepare_fastq(raw_path, intermediate_path, multithread = TRUE)
pre_primer_hit_data<- pre_primer_hit_data(primer_data, fastq_data, intermediate_path)
pre_primer_plot <- primer_hit_plot(pre_primer_hit_data, fastq_data, intermediate_path, "pre_primer_plot.pdf")
cutadapt_data <- cutadapt_tibble(fastq_data, metadata, intermediate_path)
cutadapt_run(cutadapt_path, cutadapt_data)
quality_plots<-plot_qc(cutadapt_data, intermediate_path)
filter_results <-filter_and_trim(intermediate_path, cutadapt_data, minLength=50, maxLength = 300, maxEE = 8)
post_primer_hit_data <- post_primer_hit_data(primer_data, cutadapt_data, intermediate_path)
post_primer_plot <- primer_hit_plot(post_primer_hit_data, fastq_data, intermediate_path, "post_primer_plot.pdf")
quality_plots2 <- post_trim_qc(cutadapt_data, intermediate_path)


#Command to prepare databases for downstream steps
create_ref_database(intermediate_path)
format_database(raw_path, "oomycetedb.fasta")
#If you included the ITS barcode in your analysis
#Should also test pipeline just on ITS dataset
format_database2(raw_path, "unite.fasta")
#The remaining functions will be incorporated shortly.


#TODO
#Parameters optimization
#Further clarify the inputs
#Templates for input files
#Consistency in naming
#Improve function names, and continue with documentation
#Check filter parameters-a lot of rps10 reads are being thrown out at the filter steps if you use default parameters
infer_asv_command(intermediate_path, cutadapt_data, multithread = TRUE, verbose=TRUE, nominalQ = TRUE)
merged_reads<-merge_reads_command(intermediate_path, minOverlap = 12, maxMismatch = 2, verbose=TRUE)
countOverlap(merged_reads)
raw_seqtab<-makeSeqtab(merged_reads)
asv_abund_table<-make_abund_table(raw_seqtab)
make_seqhist(raw_seqtab, 'Seqcounts_raw_abund_matrix.pdf')
make_seqhist(asv_abund_table, 'Seqcounts_postfiltered_abund_matrix.pdf')
abund_asv_rps10 <- prep_abund_table(cutadapt_data, asv_abund_table, "rps10")
abund_asv_its <- prep_abund_table(cutadapt_data, asv_abund_table, "ITS")
separate_abund_table(abund_asv_its, abund_asv_rps10, intermediate_path, asv_abund_table)
separate_abund_filepath <- file.path(intermediate_path, "Separate_abund.Rdata")
load(separate_abund_filepath)
tax_results_rps10_asv <- assign_taxonomyDada2(abund_asv_rps10, "rps10_reference_db.fa", "rps10_taxtable.Rdata", outputBootstraps = TRUE, multithread=TRUE, minBoot=0) #still neesd some work
tax_results_its_asv <- assign_taxonomyDada2(abund_asv_its, "its_short.fa", "its_taxtable.Rdata",outputBootstraps = TRUE, multithread=TRUE, minBoot=0)
its_seqs <- read_fasta(file.path(intermediate_path, "reference_databases",  'its_short.fa'))
rps10_seqs <- read_fasta(file.path(intermediate_path, "reference_databases", 'rps10_reference_db.fa'))
rps10_pids_asv <- get_pids(tax_results_rps10_asv, "rps10_reference_db.fa")
its_pids_asv <- get_pids(tax_results_its_asv, "its_short.fa")
#These functions will not work if minBoot is not 0
tax_results_rps10_asv_pid <- add_pid_to_tax(tax_results_rps10_asv, rps10_pids_asv)
tax_results_its_asv_pid <- add_pid_to_tax(tax_results_its_asv, its_pids_asv)
seq_tax_asv <- c(assignTax_as_char(tax_results_rps10_asv_pid), assignTax_as_char(tax_results_its_asv))
formatted_abund_asv<-format_abund_matrix(asv_abund_table, seq_tax_asv)
dada2_readcounts_multi_sample(asv_abund_table) #check

#For just rps10 barcode
infer_asv_command(intermediate_path, cutadapt_data, multithread = TRUE, verbose=TRUE, nominalQ = TRUE)
merged_reads<-merge_reads_command(intermediate_path, minOverlap = 12, maxMismatch = 2, verbose=TRUE)
countOverlap(merged_reads)
raw_seqtab<-makeSeqtab(merged_reads)
asv_abund_table<-make_abund_table(raw_seqtab)
make_seqhist(raw_seqtab, 'Seqcounts_raw_abund_matrix.pdf')
make_seqhist(asv_abund_table, 'Seqcounts_postfiltered_abund_matrix.pdf')
abund_asv_rps10 <- prep_abund_table(cutadapt_data, asv_abund_table, "rps10")
tax_results_rps10_asv <- assign_taxonomyDada2(abund_asv_rps10, "rps10_reference_db.fa", "rps10_taxtable.Rdata", outputBootstraps = TRUE, multithread=TRUE, minBoot=0)
rps10_pids_asv <- get_pids(tax_results_rps10_asv, "rps10_reference_db.fa")
rps10_seqs <- read_fasta(file.path(intermediate_path, "reference_databases", 'rps10_reference_db.fa'))
rps10_pids_asv <- get_pids(tax_results_rps10_asv, "rps10_reference_db.fa")
tax_results_rps10_asv_pid <- add_pid_to_tax(tax_results_rps10_asv, rps10_pids_asv)
seq_tax_asv <- assignTax_as_char(tax_results_rps10_asv_pid)
formatted_abund_asv<-format_abund_matrix(asv_abund_table, seq_tax_asv)
dada2_readcounts_multi_sample(asv_abund_table)#needs checkingung's functions


#what to do about inheritting params? How to tidy this up for greater useability? 

