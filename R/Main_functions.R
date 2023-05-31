#Prepare reads. A wrapper function to prepare reads for trimming using Cutadapt. Counts of primers on reads will be output.
#' Main command prepare reads for primer trimming
#' 
#' @param primer_path The primer data tibble created in orient_primers function
#' @param metadata_path A path to a metadata containing the concatenated metadata and primer data
#' @param fastq_path path to a directory containing FASTQ reads
#' @inheritParams read_prefilt_fastq
#' @inheritParams read_prefilt_fastq
#' @return A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @export prepare_reads
#' @examples data_tables<-prepare_reads(directory_path, primer_path, metadata_path, fastq_path, directory_path, maxN=0, multithread=TRUE)

prepare_reads <- function(directory_path, primer_path, metadata_path, fastq_path, maxN=0, multithread=FALSE){
  primer_data <- orient_primers(primer_path)
  metadata <- prepare_metadata_table(metadata_path, primer_data)
  fastq_data <- read_prefilt_fastq(directory_path, maxN = maxN, multithread = multithread)
  pre_primer_hit_data<- get_pre_primer_hits(primer_data, fastq_data, directory_path)
  pre_primer_plot <- make_primer_hit_plot(pre_primer_hit_data, fastq_data, directory_path,"pre_primer_plot.pdf")
  cutadapt_data <- make_cutadapt_tibble(fastq_data, metadata, directory_path)
  data_tables <- list(cutadapt_data=cutadapt_data, primer_data=primer_data, fastq_data=fastq_data, metadata=metadata)
  return(data_tables)
}

#Trim primers
#' Main command to trim primers based on DADA2 functions
#' 
#' @inheritParams plot_qc
#' @inheritParams filter_and_trim
#' @inheritParams post_trim_qc
#' @inheritParams run_cutadapt
#' @param data_tables A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param directory_path A path to the intermediate folder and directory
#' @param cutadapt_path A path to the cutadapt program
#' @return
#' @export cut_trim
#' @examples
cut_trim <- function(data_tables,directory_path,cutadapt_path,
                     maxEE = Inf, truncQ = 2, minLen = 20, maxLen = Inf,
                     truncLen = 0, maxN = 0, minQ=0, trimLeft=0, trimRight=0,rm.phix=TRUE,
                     multithread=FALSE, matchIDs=FALSE, verbose=FALSE,
                     qualityType="Auto", OMP=TRUE, n=1e+05,id.sep="\\s",
                     rm.lowcomplex=0, orient.fwd=NULL, id.field=NULL, minCutadaptlength=50){
  run_cutadapt(cutadapt_path, data_tables$cutadapt_data, minCutadaptlength=minCutadaptlength)
  quality_plots<-plot_qc(data_tables$cutadapt_data, directory_path)
  filter_results <- filter_and_trim(directory_path, data_tables$cutadapt_data,
                                    maxEE = maxEE, truncQ = truncQ, minLen = minLen,
                                    maxLen = maxLen, multithread = multithread,
                                    matchIDs=matchIDs, verbose = verbose, qualityType = qualityType,
                                    OMP = OMP, n=n, id.sep = id.sep, trimLeft=trimLeft, trimRight=trimRight, rm.lowcomplex = rm.lowcomplex,
                                    orient.fwd = orient.fwd, id.field = id.field)
  post_primer_hit_data <- get_post_primer_hits(data_tables$primer_data, data_tables$cutadapt_data,
                                               directory_path)
  post_primer_plot <- make_primer_hit_plot(post_primer_hit_data, data_tables$fastq_data,
                                           directory_path, "post_primer_plot.pdf")
  quality_plots2 <- plot_post_trim_qc(data_tables$cutadapt_data, directory_path)
}

#' Make an amplified sequence variants abundance matrix with data processed through preceding steps
#'
#' @param data_tables A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param rawSeqTab_fileName A filename as which the raw sequence table will be saved
#' @param abundMatrix_fileName A filenmae as which the abundance matrix will be saved
#' @param directory_path
#' @inheritParams infer_asv_command
#' @inheritParams merge_reads_command
#' @return asv dataframe
#' @export make_asv_abund_matrix
#' @examples
make_asv_abund_matrix <- function(data_tables, directory_path, cutadapt_data,
                                  multithread=FALSE,nbases = 1e+08,
                                  errorEstimationFunction = loessErrfun, randomize=FALSE,
                                  MAX_CONSIST=10, OMEGA_C=0, qualityType="Auto", nominalQ = FALSE,
                                  obs=TRUE, err_zout=TRUE, err_in=FALSE, pool=FALSE, selfConsist=FALSE,
                                  verbose=FALSE,  minOverlap=12, maxMismatch=0, returnRejects=FALSE,
                                  justConcatenate=FALSE, trimOverhang=FALSE, orderBy="abundance",
                                  method="consensus", min_asv_length=50)
{
  infer_asv_command(directory_path, data_tables$cutadapt_data,
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
  merged_reads<-merge_reads_command(directory_path, minOverlap = minOverlap,
                                    maxMismatch = maxMismatch,
                                    returnRejects = returnRejects,
                                    justConcatenate = justConcatenate,
                                    verbose = verbose)
  countOverlap(merged_reads)
  raw_seqtab<-makeSeqtab(merged_reads, orderBy=orderBy)
  asv_abund_matrix<-make_abund_matrix(raw_seqtab, min_asv_length=min_asv_length)
  make_seqhist(asv_abund_matrix)
  return(asv_abund_matrix)
}

# Hung's Addition: process_single_barcode
#' Process the information from an ASV abundance matrix to run DADA2 for single barcode
#'
#' @param data_tables A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#' @inheritParams assign_taxonomyDada2
#' @return
#' @examples
#add params
#name the ref db by barcode name
process_single_barcode <- function(data_tables, asv_abund_matrix, tryRC=FALSE, verbose=FALSE, multithread=FALSE, barcode="rps10")
{
  abund_asv_single <- prep_abund_matrix(data_tables$cutadapt_data, asv_abund_matrix, barcode)
  refdb=paste0(barcode, "_reference_db.fa")
  taxmat=paste0(barcode,"_taxmatrix.Rdata")
  tax_results_single_asv <- assign_taxonomyDada2(abund_asv_single, refdb, taxmat,
                                                 tryRC=tryRC, verbose=verbose, multithread=multithread)
  single_pids_asv <- get_pids(tax_results_single_asv, refdb)
  tax_results_single_asv_pid <- add_pid_to_tax(tax_results_single_asv, single_pids_asv)
  seq_tax_asv <- assignTax_as_char(tax_results_single_asv_pid)
  formatted_abund_asv<-format_abund_matrix(asv_abund_matrix, seq_tax_asv)
  get_read_counts(asv_abund_matrix)
}

# Hung's addition: process_rps10_ITS_barcode
#For rps10 and its barcodes
#' Main trim command to run core DADA2 functions for 2 barcodes
#'
#' @param data_tables list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#' @inheritParams assign_taxonomyDada2
#' @return
#' @examples
process_pooled_barcode <- function(data_tables, asv_abund_matrix, tryRC=FALSE, verbose=FALSE, multithread=FALSE, barcode1="rps10", barcode2="its")
{
  abund_asv_barcode1 <- prep_abund_matrix(data_tables$cutadapt_data, asv_abund_matrix, barcode1)
  abund_asv_barcode2 <- prep_abund_matrix(data_tables$cutadapt_data, asv_abund_matrix, barcode2)
  separate_abund_matrix(abund_asv_barcode2, abund_asv_barcode1, directory_path, asv_abund_matrix)
  separate_abund_filepath <- file.path(directory_path, "Separate_abund.Rdata")
  load(separate_abund_filepath)
  refdb1=paste0(barcode1, "_reference_db.fa")
  taxmat1=paste0(barcode1,"_taxmatrix.Rdata")
  refdb2=paste0(barcode2, "_reference_db.fa")
  taxmat2=paste0(barcode2,"_taxmatrix.Rdata")
  tax_results_barcode1_asv <- assign_taxonomyDada2(abund_asv_barcode1, refdb1, taxmat1,
                                                   tryRC=FALSE, verbose=FALSE, multithread=FALSE)
  tax_results_barcode2_asv <- assign_taxonomyDada2(abund_asv_barcode2, refdb2, taxmat2,
                                                   tryRC=FALSE, verbose=FALSE, multithread=FALSE)
  barcode1_pids_asv <- get_pids(tax_results_barcode1_asv, refdb1)
  barcode2_pids_asv <- get_pids(tax_results_barcode2_asv, refdb2)
  tax_results_barcode1_asv_pid <- add_pid_to_tax(tax_results_barcode1_asv, barcode1_pids_asv)
  tax_results_barcode2_asv_pid <- add_pid_to_tax(tax_results_barcode2_asv, barcode2_pids_asv)
  seq_tax_asv <- c(assignTax_as_char(tax_results_barcode1_asv_pid), assignTax_as_char(tax_results_barcode2_asv_pid))
  formatted_abund_asv<-format_abund_matrix(asv_abund_matrix, seq_tax_asv)
  get_read_counts(asv_abund_matrix)
}


#' Assign rps10 and ITS taxonomy
#'
#' @param data_tables
#' @param barcode
#' @param asv_abund_matrix
#' @param tryRC
#' @param verbose
#' @param multithread
#' @param database_rps10
#' @param database_its
#' @inheritParams assign_taxonomyDada2
#' @return
#' @export assignTax
#' @examples

#generalize for different db
assignTax <- function(directory_path, data_tables, asv_abund_matrix, tryRC=FALSE, verbose=FALSE, multithread=FALSE, barcode="rps10", database_rps10="oomycetedb.fasta", database_its="sh_general_release_dynamic_22.08.2016.fasta") {
  if(barcode=="rps10") {
    format_database_rps10(directory_path, database_rps10)
    summary_table<-process_single_barcode(data_tables, asv_abund_matrix, multithread = multithread, barcode="rps10")
  } else if(barcode=="its") {
    format_database_its(directory_path, database_its)
    summary_table<-process_single_barcode(data_tables, asv_abund_matrix, multithread = multithread, barcode="its")
  } else if(barcode=="rps10_its") {
    format_database_rps10(directory_path, database_rps10)
    format_database_its(directory_path, database_its)
    summary_table<-process_pooled_barcode(data_tables, asv_abund_matrix, multithread = multithread, barcode1="rps10", barcode2="its")
  } else {
    print("Barcodes note recognized")
  }
}

#' Filter ASV abundance matrix and convert to taxmap object
#'
#' @param asv_abund_matrix
#' @param min_read_depth
#' @param minimum_bootstrap
#' @param pid_species
#' @param pid_genus
#' @param pid_family
#' @return
#' @export asvmatrix_to_taxmap
#' @examples
asvmatrix_to_taxmap <- function(min_read_depth=0, minimum_bootstrap=50, pid_species=0, pid_genus=0, pid_family=0){
  abund_matrix <- read_csv(file.path(directory_path,'final_asv_abundance_matrix.csv'))
  is_low_abund <- rowSums(abund_matrix[, data_tables$metadata$sample_nameBarcode]) < min_read_depth
  abundance <- filter(abund_matrix, ! is_low_abund)
  pid_cutoffs <- list(species = pid_species, genus = pid_genus, family = pid_family)
  abundance$dada2_tax <- map_chr(strsplit(abundance$dada2_tax, ';'), function(x) {
    map_chr(strsplit(x, '--'), function(parts) {
      if (parts[[3]] != "ASV" && as.numeric(parts[[2]]) <= minimum_bootstrap) {
        parts[[1]] <- "Unknown"
      }
      return(paste(parts, collapse = '--'))
    }) %>% paste(collapse = ';')
  })
  obj_dada <- metacoder::parse_tax_data(abundance, class_cols = 'dada2_tax', class_sep = ';', include_tax_data = TRUE,
                                        class_regex = '^(.+)--(.+)--(.+)$',
                                        class_key = c(taxon = 'taxon_name', boot = 'info', rank = 'taxon_rank'))
  names(obj_dada$data) <- c('abund', 'score') #do I need this?
  return(obj_dada)
}

#' Convert taxmap object to Phyloseq object (metacoder wrapper function)
#'
#' @param obj_data
#' @return A phyloseq object
#' @export taxmap_to_phyloseq
#' @examples

taxmap_to_phyloseq <- function(obj_dada) {
  obj_dada$data$otu_table=obj_dada$data$abund[,-2:-4]
  obj_dada$data$otu_table$otu_id=paste0('ASV', 1:nrow(obj_dada$data$otu_table))
  obj_dada$data$sample_data=data_tables$metadata
  
  phylo_obj<-metacoder::as_phyloseq(
    obj_dada,
    sample_data = obj_dada$data$sample_data,
    sample_id_col = "sample_nameBarcode",
  )
}
