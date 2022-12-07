#Trim primers
#' Main command to trim primers based on DADA2 functions
#' @inheritParams plot_qc
#' @inheritParams filter_and_trim
#' @inheritParams post_trim_qc
#' @inheritParams run_cutadapt
#' @param returnList A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param intermediate_path A path to the intermediate folder and directory
#' @param cutadapt_path A path to the cutadapt program
#'
#' @return
#' @export
#'
#' @examples
#'
#'

cut_trim <- function(returnList,intermediate_path,cutadapt_path,
                     maxEE = Inf, truncQ = 2, minLength = 20, maxLength = Inf,
                     truncLen = 0, maxN = 0, minQ=0, rm.phix=TRUE,
                     multithread=FALSE, matchIDs=FALSE, verbose=FALSE,
                     qualityType="Auto", OMP=TRUE, n=1e+05,id.sep="\\s",
                     rm.lowcomplex=0, orient.fwd=NULL, id.field=NULL, min_length=25){
  run_cutadapt(cutadapt_path, returnList$cutadapt_data, min_length=min_length)
  quality_plots<-plot_qc(returnList$cutadapt_data, intermediate_path)
  filter_results <-filter_and_trim(intermediate_path, returnList$cutadapt_data,
                                   maxEE = maxEE, truncQ = truncQ, minLength = minLength,
                                   maxLength = maxLength, multithread = multithread,
                                   matchIDs=matchIDs, verbose = verbose, qualityType = qualityType,
                                   OMP = OMP, n=n, id.sep = id.sep, rm.lowcomplex = rm.lowcomplex,
                                   orient.fwd = orient.fwd, id.field = id.field)
  post_primer_hit_data <- get_post_primer_hits(returnList$primer_data, returnList$cutadapt_data,
                                               intermediate_path)
  post_primer_plot <- make_primer_hit_plot(post_primer_hit_data, returnList$fastq_data,
                                      intermediate_path, "post_primer_plot.pdf")
  quality_plots2 <- plot_post_trim_qc(returnList$cutadapt_data, intermediate_path)
}

# Hung's Addition: make_asvAbund_table
#' Make an amplified sequence variants abundance matrix with data processed through preceding steps
#'
#'
#' @param returnList A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param rawSeqTab_fileName A filename as which the raw sequence table will be saved
#' @param abundMatrix_fileName A filenmae as which the abundance matrix will be saved
#' @inheritParams infer_asv_command 
#' @inheritParams merge_reads_command 
#' 
#'
#' @return asv_df
#' @export
#'
#' @examples asv_abund_matrix <- make_asvAbund_matrix(returnList,intermediate_path, returnList$cutadapt_data)
make_asv_abund_matrix <- function(returnList, intermediate_path, cutadapt_data,
                                 multithread=FALSE,nbases = 1e+08,
                                 errorEstimationFunction = loessErrfun, randomize=FALSE,
                                 MAX_CONSIST=10, OMEGA_C=0, qualityType="Auto", nominalQ = FALSE, 
                                 obs=TRUE, err_out=TRUE, err_in=FALSE, pool=FALSE, selfConsist=FALSE, 
                                 verbose=FALSE,  minOverlap=12, maxMismatch=0, returnRejects=FALSE, 
                                 justConcatenate=FALSE, trimOverhang=FALSE, orderBy="abundance",
                                 method="consensus", min_asv_length=50)
{
  infer_asv_command(intermediate_path, returnList$cutadapt_data, 
                    multithread = multithread, 
                    nbases=nbases,
                    errorEstimationFunction= errorEstimationFunction, 
                    randomize=randomize, 
                    MAX_CONSIST=MAX_CONSIST, 
                    OMEGA_C= OMEGA_C, 
                    qualityType=qualityType, 
                    nominalQ = nominalQ, 
                    obs=obs, 
                    err_out=err_out, 
                    err_in=err_in,
                    pool=pool, 
                    selfConsist=selfConsist, 
                    verbose=verbose)
  merged_reads<-merge_reads_command(intermediate_path, minOverlap = minOverlap,
                                    maxMismatch = maxMismatch,
                                    returnRejects = returnRejects,
                                    justConcatenate = justConcatenate,
                                    verbose = verbose)
  countOverlap(merged_reads)
  raw_seqtab<-makeSeqtab(merged_reads, orderBy=orderBy)
  asv_abund_matrix<-make_abund_table(raw_seqtab, min_asv_length=min_asv_length)
  make_seqhist(asv_abund_matrix)
  return(asv_abund_matrix)
}

# Hung's Addition: process_rps10_barcode
#' Process the information from an ASV abundance matrix to run DADA2 for rps10 barcode
#'
#'
#' @param returnList A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#' @inheritParams assign_taxonomyDada2 
#' @return 
#' @export
#'
#' @examples
#add params
process_rps10_barcode <- function(returnList, asv_abund_matrix, tryRC=FALSE, verbose=FALSE, multithread=FALSE)
{
  abund_asv_rps10 <- prep_abund_table(returnList$cutadapt_data, asv_abund_matrix, "rps10")
  tax_results_rps10_asv <- assign_taxonomyDada2(abund_asv_rps10, "rps10_reference_db.fa", "rps10_taxtable.Rdata",
                                                tryRC=tryRC, verbose=verbose, multithread=multithread)
  rps10_pids_asv <- get_pids(tax_results_rps10_asv, "rps10_reference_db.fa")
  tax_results_rps10_asv_pid <- add_pid_to_tax(tax_results_rps10_asv, rps10_pids_asv)
  seq_tax_asv <- assignTax_as_char(tax_results_rps10_asv_pid)
  formatted_abund_asv<-format_abund_matrix(asv_abund_matrix, seq_tax_asv)
  get_read_counts(asv_abund_matrix)
}

# Hung's addition: process_rps10_ITS_barcode
#For rps10 and its barcodes
#' Main trim command to run core DADA2 functions for 2 barcodes
#'
#' @param returnList list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#' @inheritParams assign_taxonomyDada2 
#' 
#' @return
#' @export
#'
#' @examples
process_rps10_ITS_barcode <- function(returnList, asv_abund_matrix, tryRC=FALSE, verbose=FALSE, multithread=FALSE)
{
  abund_asv_rps10 <- prep_abund_table(returnList$cutadapt_data, asv_abund_matrix, "rps10")
  abund_asv_its <- prep_abund_table(returnList$cutadapt_data, asv_abund_matrix, "ITS")
  separate_abund_table(abund_asv_its, abund_asv_rps10, intermediate_path, asv_abund_matrix)
  separate_abund_filepath <- file.path(intermediate_path, "Separate_abund.Rdata")
  load(separate_abund_filepath)
  tax_results_rps10_asv <- assign_taxonomyDada2(abund_asv_rps10, "rps10_reference_db.fa", "rps10_taxtable.Rdata",
                                                tryRC=FALSE, verbose=FALSE, multithread=FALSE) #FIX
  tax_results_its_asv <- assign_taxonomyDada2(abund_asv_its, "its_short.fasta", "its_taxtable.Rdata",
                                              tryRC=FALSE, verbose=FALSE, multithread=FALSE)
  rps10_pids_asv <- get_pids(tax_results_rps10_asv, "rps10_reference_db.fa")
  its_pids_asv <- get_pids(tax_results_its_asv, "its_short.fasta")
  tax_results_rps10_asv_pid <- add_pid_to_tax(tax_results_rps10_asv, rps10_pids_asv)
  tax_results_its_asv_pid <- add_pid_to_tax(tax_results_its_asv, its_pids_asv)
  seq_tax_asv <- c(assignTax_as_char(tax_results_rps10_asv_pid), assignTax_as_char(tax_results_its_asv_pid))
  formatted_abund_asv<-format_abund_matrix(asv_abund_matrix, seq_tax_asv)
  get_read_counts(asv_abund_matrix)
}
