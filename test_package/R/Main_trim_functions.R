#' A function for calling read_fastq, matching_order_primer_check, and remove_ns functions. This will process and edit the FASTQ and make them ready for the trimming of primers with Cutadapt.
#' @inheritParams dada2::filterAndTrim
#' @param raw_path A path to a directory that contains raw data
#' @param intermediate_path A path to the intermediate folder and directory
#' @param metadata A metadata containing the concatenated metadata and primer data

#'
#' @return
#' @export #add params
#'
#' @examples
prepare_fastq <- function(raw_path,intermediate_path, maxN= 0, multithread = FALSE){
  fastq_data <- read_fastq(raw_path)
  matching_order_primer_check(fastq_data)
  fastq_data <- remove_ns(fastq_data, intermediate_path, maxN, multithread=multithread)

  return(fastq_data)
}
#' Runs all of the functions associated with running cutadapt. Users should not have to touch these functions.
#'
#' @param directory_path A path to a directory containing reads and metadata/primer files
#' @param primer_path The primer data tibble created in prepare_primers function
#' @param metadata_path A path to a metadata containing the concatenated metadata and primer data
#' @param fastq_path A path to a directory containing FASTQ reads
#' @param intermediate_path A path to the intermediate folder and directory
#' @param cutadapt_path A path to the cutadapt program
#'
#' @return Intermediate_data folder with trimmed and filtered reads for each sample.
#' @export
#'
#' @examples
main_cutadapt_function <- function(directory_path, primer_path, metadata_path, fastq_path,intermediate_path, cutadapt_path){
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
}

#For more than 1 sample, rps10 barcode
#' Main trim command to run core DADA2 functions for rps10 barcode
#'
#' @param intermediate_path A path to the intermediate folder and directory
#' @param cutadapt_data Intermediate_data folder with trimmed and filtered reads for each sample
#'
#' @return
#' @export
#'
#' @examples
rps10_barcode_function<- function(intermediate_path, cutadapt_data){
  infer_asv_command(intermediate_path, cutadapt_data)
  merged_reads<-merge_reads_command(intermediate_path)
  countOverlap(merged_reads)
  raw_seqtab<-makeSeqtab(merged_reads)
  asv_abund_table<-make_abund_table(raw_seqtab)
  make_seqhist(raw_seqtab, 'Seqcounts_raw_abund_matrix.pdf')
  make_seqhist(asv_abund_table, 'Seqcounts_postfiltered_abund_matrix.pdf')
  abund_asv_rps10 <- prep_abund_table(cutadapt_data, asv_abund_table, "rps10")
  tax_results_rps10_asv <- assign_taxonomyDada2(abund_asv_rps10, "rps10_reference_db.fa", "rps10_taxtable.Rdata")
  rps10_pids_asv <- get_pids(tax_results_rps10_asv, "rps10_reference_db.fa")
  tax_results_rps10_asv_pid <- add_pid_to_tax(tax_results_rps10_asv, rps10_pids_asv)
  seq_tax_asv <- assignTax_as_char(tax_results_rps10_asv_pid)
  formatted_abund_asv<-format_abund_matrix(asv_abund_table, seq_tax_asv)
  dada2_readcounts_multi_sample(formatted_abund_asv) #needs checking
}


#For rps10 and its barcodes
#' Main trim command to run core DADA2 functions for 2 barcodes
#'
#' @param intermediate_path A path to the intermediate folder and directory
#' @param cutadapt_data Intermediate_data folder with trimmed and filtered reads for each sample
#'
#' @return
#' @export
#'
#' @examples
ITS_rps10_barcode_function<- function(intermediate_path, cutadapt_data){
  infer_asv_command(intermediate_path, cutadapt_data)
  merged_reads<-merge_reads_command(intermediate_path)
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
  tax_results_rps10_asv <- assign_taxonomyDada2(abund_asv_rps10, "rps10_reference_db.fa", "rps10_taxtable.Rdata") #still neesd some work
  tax_results_its_asv <- assign_taxonomyDada2(abund_asv_its, "its_short.fasta", "its_taxtable.Rdata")
  its_seqs <- read_fasta(file.path(intermediate_path, "reference_databases",  'its_short.fasta'))
  rps10_seqs <- read_fasta(file.path(intermediate_path, "reference_databases", 'rps10_reference_db.fa'))
  rps10_pids_asv <- get_pids(tax_results_rps10_asv, "rps10_reference_db.fa")
  its_pids_asv <- get_pids(tax_results_its_asv, "its_short.fasta")
  tax_results_rps10_asv_pid <- add_pid_to_tax(tax_results_rps10_asv, rps10_pids_asv)
  tax_results_its_asv_pid <- add_pid_to_tax(tax_results_its_asv, its_pids_asv)
  seq_tax_asv <- c(assignTax_as_char(tax_results_rps10_asv_pid), assignTax_as_char(tax_results_its_asv))
  formatted_abund_asv<-format_abund_matrix(asv_abund_table, seq_tax_asv)
  dada2_readcounts_multi_sample(formatted_abund_asv) #check
}


#Hung's additions
#Prepare primers
#' Main command prepare reads for primer trimming
#' @param directory_path A path to a directory containing reads and metadata/primer files
#' @param primer_path The primer data tibble created in prepare_primers function
#' @param metadata_path A path to a metadata containing the concatenated metadata and primer data
#' @param fastq_path path to a directory containing FASTQ reads
#' @param fastq_pathA path to a directory containing FASTQ reads
#' @param intermediate_path A path to the intermediate folder and directory
#' @inheritParams prepare_fastq
#'
#' @return A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @export
#'
#' @examples

prepare <- function(directory_path, primer_path, metadata_path, fastq_path,intermediate_path, maxN=0, multithread=FALSE){
  intermediate_path <- create_intermediate(directory_path)
  prep_tables <- file.path(intermediate_path, "Prep_tables.Rdata")
  primer_data <- prepare_primers(primer_path)
  metadata <- prepare_metadata(metadata_path, primer_data)
  fastq_data <- prepare_fastq(raw_path, intermediate_path, maxN = maxN, multithread = multithread)
  pre_primer_hit_data<- pre_primer_hit_data(primer_data, fastq_data,
                                            intermediate_path)
  pre_primer_plot <- primer_hit_plot(pre_primer_hit_data, fastq_data,
                                     intermediate_path, "pre_primer_plot.pdf")
  cutadapt_data <- cutadapt_tibble(fastq_data, metadata, intermediate_path)
  returnList <- list(cutadapt_data=cutadapt_data, primer_data=primer_data, fastq_data=fastq_data, metadata=metadata)
  return(returnList)

}

#Hung's additions
#Trim primers
#' Main command to trim primers based on DADA2 functions
#' @inheritParams plot_qc
#' @inheritParams filter_and_trim
#' @inheritParams post_trim_qc
#' @param returnList A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param intermediate_path A path to the intermediate folder and directory
#' @param cutadapt_path A path to the cutadapt program
#'
#' @return
#' @export
#'
#' @examples
#'
cut_trim <- function(returnList,intermediate_path,cutadapt_path,
                     maxEE = Inf, truncQ = 2, minLength = 20, maxLength = Inf,
                     truncLen = 0, maxN = 0, minQ=0, rm.phix=TRUE,
                     multithread=FALSE, verbose=FALSE,
                     n=1e+05){
  cutadapt_run(cutadapt_path, returnList$cutadapt_data)
  quality_plots<-plot_qc(returnList$cutadapt_data, intermediate_path)
  filter_results <-filter_and_trim(intermediate_path, returnList$cutadapt_data,
                                   maxEE = maxEE, truncQ = truncQ, minLength = minLength,
                                   maxLength = maxLength, truncLen=truncLen, maxN=maxN, minQ=minQ, 
                                   multithread = multithread, rm.phix=rm.phix, 
                                   verbose = verbose, n=n)
  post_primer_hit_data <- post_primer_hit_data(returnList$primer_data, returnList$cutadapt_data,
                                               intermediate_path)
  post_primer_plot <- primer_hit_plot(post_primer_hit_data, returnList$fastq_data,
                                      intermediate_path, "post_primer_plot.pdf")
  quality_plots2 <- post_trim_qc(returnList$cutadapt_data, intermediate_path)
}


# Hung's Addition: make_asvAbund_table
#' Make an amplified sequence variants abundance matrix with data processed through preceding steps
#'
#'
#' @param returnList A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param rawSeqTab_fileName A filename as which the raw sequence table will be saved
#' @param abundMatrix_fileName A filenmae as which the abundance matrix will be saved
#'
#' @return asv_df
#' @export
#'
#' @examples asv_abund_matrix <- make_asvAbund_matrix(returnList)
# Clarify slight restructuring with Hung
make_asvAbund_matrix <- function(returnList, intermediate_path){
  infer_asv_command(intermediate_path, cutadapt_data=returnList$cutadapt_data)
  merged_reads<-merge_reads_command(intermediate_path)
  countOverlap(merged_reads)
  raw_seqtab<-makeSeqtab(merged_reads)
  asv_abund_matrix<-make_abund_table(raw_seqtab)
  make_seqhist(asv_abund_matrix)
  return(asv_abund_matrix)
}

# Hung's Addition: process_rps10_barcode
#' Process the information from an ASV abundance matrix to run DADA2 for rps10 barcode
#'
#'
#' @param returnList A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#'
#' @return
#' @export
#'
#' @examples
process_rps10_barcode <- function(returnList, asv_abund_matrix)
{
  abund_asv_rps10 <- prep_abund_table(returnList$cutadapt_data, asv_abund_matrix, "rps10")
  tax_results_rps10_asv <- assign_taxonomyDada2(abund_asv_rps10, "rps10_reference_db.fa", "rps10_taxtable.Rdata")
  rps10_pids_asv <- get_pids(tax_results_rps10_asv, "rps10_reference_db.fa")
  tax_results_rps10_asv_pid <- add_pid_to_tax(tax_results_rps10_asv, rps10_pids_asv)
  seq_tax_asv <- assignTax_as_char(tax_results_rps10_asv_pid)
  formatted_abund_asv<-format_abund_matrix(asv_abund_table, seq_tax_asv)
  dada2_readcounts_multi_sample(formatted_abund_asv)
}

# Hung's addition: process_rps10_ITS_barcode
#For rps10 and its barcodes
#' Main trim command to run core DADA2 functions for 2 barcodes
#'
#' @param returnListA list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#'
#' @return
#' @export
#'
#' @examples
process_rps10_ITS_barcode <- function(returnList, intermediate_path, asv_abund_matrix)
{
  abund_asv_rps10 <- prep_abund_table(cutadapt_data, asv_abund_matrix, "rps10")
  abund_asv_its <- prep_abund_table(cutadapt_data, asv_abund_matrix, "ITS")
  separate_abund_table(abund_asv_its, abund_asv_rps10, intermediate_path, asv_abund_matrix)
  separate_abund_filepath <- file.path(intermediate_path, "Separate_abund.Rdata")
  load(separate_abund_filepath)
  tax_results_rps10_asv <- assign_taxonomyDada2(abund_asv_rps10, "rps10_reference_db.fa", "rps10_taxtable.Rdata") #still neesd some work
  tax_results_its_asv <- assign_taxonomyDada2(abund_asv_its, "its_short.fasta", "its_taxtable.Rdata")
  its_seqs <- read_fasta(file.path(intermediate_path, "reference_databases",  'its_short.fasta'))
  rps10_seqs <- read_fasta(file.path(intermediate_path, "reference_databases", 'rps10_reference_db.fa'))
  rps10_pids_asv <- get_pids(tax_results_rps10_asv, "rps10_reference_db.fa")
  its_pids_asv <- get_pids(tax_results_its_asv, "its_short.fasta")
  tax_results_rps10_asv_pid <- add_pid_to_tax(tax_results_rps10_asv, rps10_pids_asv)
  tax_results_its_asv_pid <- add_pid_to_tax(tax_results_its_asv, its_pids_asv)
  seq_tax_asv <- c(assignTax_as_char(tax_results_rps10_asv_pid), assignTax_as_char(tax_results_its_asv))
  formatted_abund_asv<-format_abund_matrix(asv_abund_matrix, seq_tax_asv)
  dada2_readcounts_multi_sample(formatted_abund_asv)
}
#Hung's additions
#Trim primers
#' Main command to trim primers based on DADA2 functions
#' @inheritParams plot_qc
#' @inheritParams filter_and_trim
#' @inheritParams post_trim_qc
#' @param returnList A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param intermediate_path A path to the intermediate folder and directory
#' @param cutadapt_path A path to the cutadapt program
#'
#' @return
#' @export
#'
#' @examples
#'
cut_trim <- function(returnList,intermediate_path,cutadapt_path,
                     maxEE = Inf, truncQ = 2, minLength = 20, maxLength = Inf,
                     truncLen = 0, maxN = 0, minQ=0, rm.phix=TRUE,
                     multithread=FALSE, matchIDs=FALSE, verbose=FALSE,
                     qualityType="Auto", OMP=TRUE, n=1e+05,id.sep="\\s",
                     rm.lowcomplex=0, orient.fwd=NULL, id.field=NULL){
  cutadapt_run(cutadapt_path, returnList$cutadapt_data)
  quality_plots<-plot_qc(returnList$cutadapt_data, intermediate_path)
  filter_results <-filter_and_trim(intermediate_path, returnList$cutadapt_data,
                                   maxEE = maxEE, truncQ = truncQ, minLength = minLength,
                                   maxLength = maxLength, multithread = multithread,
                                   matchIDs=matchIDs, verbose = verbose, qualityType = qualityType,
                                   OMP = OMP, n=n, id.sep = id.sep, rm.lowcomplex = rm.lowcomplex,
                                   orient.fwd = orient.fwd, id.field = id.field)
  post_primer_hit_data <- post_primer_hit_data(returnList$primer_data, returnList$cutadapt_data,
                                               intermediate_path)
  post_primer_plot <- primer_hit_plot(post_primer_hit_data, returnList$fastq_data,
                                      intermediate_path, "post_primer_plot.pdf")
  quality_plots2 <- post_trim_qc(returnList$cutadapt_data, intermediate_path)
}


# Hung's Addition: make_asvAbund_table
#' Make an amplified sequence variants abundance matrix with data processed through preceding steps
#'
#'
#' @param returnList A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param rawSeqTab_fileName A filename as which the raw sequence table will be saved
#' @param abundMatrix_fileName A filenmae as which the abundance matrix will be saved
#'
#' @return asv_df
#' @export
#'
#' @examples asv_abund_matrix <- make_asvAbund_matrix(returnList,intermediate_path, returnList$cutadapt_data)
make_asvAbund_matrix <- function(returnList, intermediate_path, cutadapt_data,
                                rawSeqTab_fileName = 'Seqcounts_raw_abund_matrix.pdf',
                                abundMatrix_fileName = 'Seqcounts_postfiltered_abund_matrix.pdf')
{
  infer_asv_command(intermediate_path, returnList$cutadapt_data)
  merged_reads<-merge_reads_command(intermediate_path)
  countOverlap(merged_reads)
  raw_seqtab<-makeSeqtab(merged_reads)
  asv_abund_matrix<-make_abund_table(raw_seqtab)
  make_seqhist(asv_abund_matrix)
  return(asv_abund_matrix)
}


# Hung's Addition: process_rps10_barcode
#' Process the information from an ASV abundance matrix to run DADA2 for rps10 barcode
#'
#'
#' @param returnList A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#'
#' @return
#' @export
#'
#' @examples
process_rps10_barcode <- function(returnList, asv_abund_matrix)
{
  abund_asv_rps10 <- prep_abund_table(ReturnList$cutadapt_data, asv_abund_matrix, "rps10")
  tax_results_rps10_asv <- assign_taxonomyDada2(abund_asv_rps10, "rps10_reference_db.fa", "rps10_taxtable.Rdata")
  rps10_pids_asv <- get_pids(tax_results_rps10_asv, "rps10_reference_db.fa")
  tax_results_rps10_asv_pid <- add_pid_to_tax(tax_results_rps10_asv, rps10_pids_asv)
  seq_tax_asv <- assignTax_as_char(tax_results_rps10_asv_pid)
  formatted_abund_asv<-format_abund_matrix(asv_abund_matrix, seq_tax_asv)
  dada2_readcounts_multi_sample(formatted_abund_asv)
}

# Hung's addition: process_rps10_ITS_barcode
#For rps10 and its barcodes
#' Main trim command to run core DADA2 functions for 2 barcodes
#'
#' @param returnListA list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#'
#' @return
#' @export
#'
#' @examples
process_rps10_ITS_barcode <- function(returnList, intermediate_path, asv_abund_matrix)
{
  abund_asv_rps10 <- prep_abund_table(returnList$cutadapt_data, asv_abund_matrix, "rps10")
  abund_asv_its <- prep_abund_table(returnList$cutadapt_data, asv_abund_matrix, "ITS")
  separate_abund_table(abund_asv_its, abund_asv_rps10, intermediate_path, asv_abund_matrix)
  separate_abund_filepath <- file.path(intermediate_path, "Separate_abund.Rdata")
  load(separate_abund_filepath)
  tax_results_rps10_asv <- assign_taxonomyDada2(abund_asv_rps10, "rps10_reference_db.fa", "rps10_taxtable.Rdata") #still neesd some work
  tax_results_its_asv <- assign_taxonomyDada2(abund_asv_its, "its_short.fasta", "its_taxtable.Rdata")
  its_seqs <- read_fasta(file.path(intermediate_path, "reference_databases",  'its_short.fasta'))
  rps10_seqs <- read_fasta(file.path(intermediate_path, "reference_databases", 'rps10_reference_db.fa'))
  rps10_pids_asv <- get_pids(tax_results_rps10_asv, "rps10_reference_db.fa")
  its_pids_asv <- get_pids(tax_results_its_asv, "its_short.fasta")
  tax_results_rps10_asv_pid <- add_pid_to_tax(tax_results_rps10_asv, rps10_pids_asv)
  tax_results_its_asv_pid <- add_pid_to_tax(tax_results_its_asv, its_pids_asv)
  seq_tax_asv <- c(assignTax_as_char(tax_results_rps10_asv_pid), assignTax_as_char(tax_results_its_asv))
  formatted_abund_asv<-format_abund_matrix(asv_abund_matrix, seq_tax_asv)
  dada2_readcounts_multi_sample(formatted_abund_asv)
}

#cutadapt_run(cutadapt_path, cutadapt_data)
#quality_plots<-plot_qc(cutadapt_data, intermediate_path)
#filter_results <-filter_and_trim(intermediate_path, cutadapt_data, minLength=50, maxLength = 300, maxEE = 8)
#post_primer_hit_data <- post_primer_hit_data(primer_data, cutadapt_data, intermediate_path)
#post_primer_plot <- primer_hit_plot(post_primer_hit_data, fastq_data, intermediate_path, "post_primer_plot.pdf")
#quality_plots2 <- post_trim_qc(cutadapt_data, intermediate_path)
