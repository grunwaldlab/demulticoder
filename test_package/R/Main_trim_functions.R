#' A function for calling read_fastq, matching_order_primer_check, and remove_ns functions. This will process and edit the FASTQ and make them ready for the trimming of primers with Cutadapt.
#'
#' @param raw_path
#' @param intermediate_path
#' @param metadata
#'
#' @return
#' @export
#'
#' @examples
prepare_fastq <- function(raw_path,intermediate_path){
  fastq_data <- read_fastq(raw_path)
  matching_order_primer_check(fastq_data)
  fastq_data <- remove_ns(fastq_data, intermediate_path)

  return(fastq_data)
}
#' Runs all of the functions associated with running cutadapt. Users should not have to touch these functions.
#'
#' @param directory_path
#' @param primer_path
#' @param metadata_path
#' @param fastq_path
#' @param intermediate_path
#' @param cutadapt_path
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
  dada2_readcounts_multi_sample(asv_abund_table) #needs checking
}


#For rps10 and its barcodes
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
  dada2_readcounts_multi_sample(asv_abund_table) #check
}
