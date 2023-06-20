#' Assign rps10 and ITS taxonomy
#'
#' @param directory_path Location of read files and metadata file
#' @param data_tables specify dataframe
#' @param barcode specify which barcode you have used, 'rps10', 'its', or 'rps10_its'
#' @param asv_abund_matrix specify the ASV abundance matrix for which taxonomic assignments will be given
#' @inheritParams assign_taxonomyDada2
#' @inheritParams dada2::assignTaxonomy
#' @return Taxonomic assignments of each unique ASV sequence
#' @export assignTax
#' @examples
#' directory_path<-"~/rps10package/raw_data/rps10_ITS"
#' primer_path <-file.path(directory_path, "primer_info.csv")
#' metadata_path <-file.path(directory_path,"metadata.csv")
#' cutadapt_path<-"/opt/homebrew/bin/cutadapt"
#' data_tables <-
#' prepare_reads(
#' directory_path,
#' primer_path,
#' metadata_path,
#' maxN = 0
#' )
#' cut_trim(
#' directory_path,
#' cutadapt_path,
#' verbose = TRUE,
#' maxEE = 2,
#' truncQ = 5,
#' minLen = 200,
#' maxLen = 297,
#' minCutadaptlength = 50
#')
#' asv_abund_matrix <-
#' make_asv_abund_matrix(
#' directory_path,
#' minOverlap = 15,
#' maxMismatch = 2,
#' verbose = TRUE
#' )
#' summary <- assignTax(
#' directory_path,
#' data_tables,
#' asv_abund_matrix,
#' multithread = TRUE,
#' barcode = "rps10"
#' )
#'
assignTax <-
  function(directory_path,
           data_tables,
           asv_abund_matrix,
           tryRC = FALSE,
           verbose = FALSE,
           multithread = FALSE,
           barcode = "rps10",
           rps10_db = "oomycetedb.fasta",
           its_db = "fungidb.fasta") {
    if (barcode == "rps10") {
      format_database_rps10(directory_path, rps10_db)
      summary_table <-
        process_single_barcode(data_tables,
                               asv_abund_matrix,
                               multithread = multithread,
                               barcode = "rps10")
    } else if (barcode == "its") {
      format_database_its(directory_path, its_db)
      summary_table <-
        process_single_barcode(data_tables,
                               asv_abund_matrix,
                               multithread = multithread,
                               barcode = "its")
    } else if (barcode == "rps10_its") {
      format_database_rps10(directory_path, rps10_db)
      format_database_its(directory_path, its_db)
      summary_table <-
        process_pooled_barcode(
          data_tables,
          asv_abund_matrix,
          multithread = multithread,
          barcode1 = "rps10",
          barcode2 = "its"
        )
    } else {
      print("Barcodes not recognized")
    }
  }

#' Process the information from an ASV abundance matrix to run DADA2 for single barcode
#'
#' @param data_tables A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#' @inheritParams assign_taxonomyDada2
#' @inheritParams dada2::assignTaxonomy
#' @keywords internal
#add params
#name the ref db by barcode name
process_single_barcode <-
  function(data_tables,
           asv_abund_matrix,
           tryRC = FALSE,
           verbose = FALSE,
           multithread = FALSE,
           barcode = "rps10")
  {
    abund_asv_single <-
      prep_abund_matrix(data_tables$cutadapt_data, asv_abund_matrix, barcode)
    refdb = paste0(barcode, "_reference_db.fa")
    taxmat = paste0(barcode, "_taxmatrix.Rdata")
    tax_results_single_asv <-
      assign_taxonomyDada2(
        abund_asv_single,
        refdb,
        taxmat,
        tryRC = tryRC,
        verbose = verbose,
        multithread = multithread
      )
    single_pids_asv <- get_pids(tax_results_single_asv, refdb)
    tax_results_single_asv_pid <-
      add_pid_to_tax(tax_results_single_asv, single_pids_asv)
    seq_tax_asv <- assignTax_as_char(tax_results_single_asv_pid)
    formatted_abund_asv <-
      format_abund_matrix(asv_abund_matrix, seq_tax_asv)
    get_read_counts(asv_abund_matrix)
  }

#' Main trim command to run core DADA2 functions for 2 barcodes
#'
#' @param data_tables list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#' @inheritParams assign_taxonomyDada2
#' @keywords internal
process_pooled_barcode <-
  function(data_tables,
           asv_abund_matrix,
           tryRC = FALSE,
           verbose = FALSE,
           multithread = FALSE,
           barcode1 = "rps10",
           barcode2 = "its")
  {
    abund_asv_barcode1 <-
      prep_abund_matrix(data_tables$cutadapt_data, asv_abund_matrix, barcode1)
    abund_asv_barcode2 <-
      prep_abund_matrix(data_tables$cutadapt_data, asv_abund_matrix, barcode2)
    separate_abund_matrix(abund_asv_barcode2,
                          abund_asv_barcode1,
                          directory_path,
                          asv_abund_matrix)
    separate_abund_filepath <-
      file.path(directory_path, "Separate_abund.Rdata")
    load(separate_abund_filepath)
    refdb1 = paste0(barcode1, "_reference_db.fa")
    taxmat1 = paste0(barcode1, "_taxmatrix.Rdata")
    refdb2 = paste0(barcode2, "_reference_db.fa")
    taxmat2 = paste0(barcode2, "_taxmatrix.Rdata")
    tax_results_barcode1_asv <-
      assign_taxonomyDada2(
        abund_asv_barcode1,
        refdb1,
        taxmat1,
        tryRC = FALSE,
        verbose = FALSE,
        multithread = FALSE
      )
    tax_results_barcode2_asv <-
      assign_taxonomyDada2(
        abund_asv_barcode2,
        refdb2,
        taxmat2,
        tryRC = FALSE,
        verbose = FALSE,
        multithread = FALSE
      )
    barcode1_pids_asv <- get_pids(tax_results_barcode1_asv, refdb1)
    barcode2_pids_asv <- get_pids(tax_results_barcode2_asv, refdb2)
    tax_results_barcode1_asv_pid <-
      add_pid_to_tax(tax_results_barcode1_asv, barcode1_pids_asv)
    tax_results_barcode2_asv_pid <-
      add_pid_to_tax(tax_results_barcode2_asv, barcode2_pids_asv)
    seq_tax_asv <-
      c(
        assignTax_as_char(tax_results_barcode1_asv_pid),
        assignTax_as_char(tax_results_barcode2_asv_pid)
      )
    formatted_abund_asv <-
      format_abund_matrix(asv_abund_matrix, seq_tax_asv)
    get_read_counts(asv_abund_matrix)
  }

#' Prepare final ASV abundance matrix
#'
#' @param directory_data folder with trimmed and filtered reads for each sample
#' @param asv_abund_matrix The returned final ASV abundance matrix
#' @param locus The barcode selected in the analysis
#' @keywords internal
prep_abund_matrix <-function(cutadapt_data, asv_abund_matrix, locus){
  #rownames(asv_abund_matrix) <- sub(rownames(asv_abund_matrix), pattern = ".fastq.gz", replacement = "")
  data_tables$cutadapt_data$file_id_primer <- paste0(data_tables$cutadapt_data$sample_name, "_", data_tables$cutadapt_data$primer_name)
  asv_abund_matrix<- asv_abund_matrix[rownames(asv_abund_matrix) %in% data_tables$cutadapt_data$file_id_primer[data_tables$cutadapt_data$primer_name == locus], ]
  return(asv_abund_matrix)
}

#' Separate abundance matrices
#'
#' @param abund_asv_barcode2 The separated matrix for second barcode 
#' @param abund_asv_barcode1 The separated matrix for first barcode
#' @param directory_path A path to the intermediate folder and directory
#' @param asv_abund_matrix The non-separated ASV matrix 
#' @keywords internal
separate_abund_matrix <- function(abund_asv_barcode2, abund_asv_barcode1, directory_path, asv_abund_matrix){
  separate_abund_path <- file.path(directory_path, "Separate_abund.Rdata")
  in_both <- colSums(abund_asv_barcode2) != 0 & colSums(abund_asv_barcode1) != 0
  assign_to_barcode2 <- in_both & colSums(abund_asv_barcode1) > colSums(abund_asv_barcode2)
  assign_to_barcode1 <- in_both & colSums(abund_asv_barcode1) < colSums(abund_asv_barcode2)
  is_barcode2 <- (colSums(abund_asv_barcode2) != 0 & colSums(abund_asv_barcode1) == 0) | assign_to_barcode2
  is_barcode1 <- (colSums(abund_asv_barcode1) != 0 & colSums(abund_asv_barcode2) == 0) | assign_to_barcode1
  abund_asv_barcode2 <- abund_asv_barcode2[ , is_barcode2]
  abund_asv_barcode1 <- abund_asv_barcode1[ , is_barcode1]
  #The number of ASVs left in the two groups should sum to the total number of ASVs, since there should be no overlap.
  save(abund_asv_barcode2, abund_asv_barcode1, file = separate_abund_path)
  stopifnot(ncol(abund_asv_barcode2) + ncol(abund_asv_barcode1) == ncol(asv_abund_matrix)) #make more messaging
}
#' Assign taxonomy
#'
#' @inheritParams dada2::assignTaxonomy
#' @param asv_abund_matrix The ASV abundance matrix
#' @param ref_database The reference database used for taxonomic inference steps
#' @param taxresults_file The name of the file for saving taxonomic assignment results
#' @keywords internal
assign_taxonomyDada2<-function(asv_abund_matrix, ref_database, taxresults_file, minBoot=0, tryRC=FALSE, verbose=FALSE, multithread=TRUE){
  tax_results<- dada2::assignTaxonomy(asv_abund_matrix,
                                      refFasta = file.path(directory_path, ref_database),
                                      taxLevels = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Reference"),
                                      minBoot = minBoot,
                                      tryRC = tryRC,
                                      outputBootstraps = TRUE,
                                      multithread = multithread)
  #if statement for bacteria
  tax_matrix_path <- file.path(directory_path, taxresults_file)
  save(tax_results, file = tax_matrix_path)
  return(tax_results)
  #save these results
}

#' Align ASV sequences to reference sequences from database to get percent ID. Get percent identities.
#'
#' @param tax_results The dataframe containing taxonomic assignments
#' @keywords internal
get_pids <- function(tax_results, db) {
  db_seqs <- read_fasta(file.path(directory_path, db))
  ref_seqs <- get_ref_seq(tax_results, db_seqs)
  asv_seqs <- rownames(tax_results$tax)
  asv_ref_align_normal <- map2(asv_seqs, ref_seqs, pairwiseAlignment, type = 'global-local')
  asv_pids_normal <- map_dbl(asv_ref_align_normal, pid, type = "PID1")
  asv_ref_align_revcomp <- map2(asv_seqs, rev_comp(ref_seqs), pairwiseAlignment, type = 'global-local')
  asv_pids_revcomp <- map_dbl(asv_ref_align_revcomp, pid, type = "PID1")
  asv_ref_align <- ifelse(asv_pids_normal >= asv_pids_revcomp, asv_ref_align_normal, asv_ref_align_revcomp)
  asv_pids <- ifelse(asv_pids_normal >= asv_pids_revcomp, asv_pids_normal, asv_pids_revcomp)
  paste0(map_chr(seq_along(asv_ref_align), function(i) {
    paste0('asv: ',  asv_seqs[[i]], '\n',
           'ref: ', names(ref_seqs)[[i]], '\n',
           'pid: ', asv_pids[[i]], '\n',
           asv_ref_align[[i]]@pattern, '\n',
           asv_ref_align[[i]]@subject, '\n')
  }), collapse = '==================================================\n') %>%
    write_lines(file = file.path(directory_path, 'dada2_asv_alignments.txt'))
  return(asv_pids)
}

#' Align ASV sequences to reference sequences from database to get percent ID. STart by retrieving reference sequences.
#'
#' @param tax_results The dataframe containing taxonomic assignments
#' @param db The reference database
#' @keywords internal
get_ref_seq <- function(tax_results, db) {
  ref_i <- as.integer(str_match(tax_results$tax[, 'Reference'], '^.+_([0-9]+)$')[ ,2])
  db[ref_i]
}

#' Add PID and bootstrap values to tax result.
#'
#' @param tax_results The dataframe containing taxonomic assignments
#' @param asv_pid Percent identity information for each ASV relative to reference database sequence
#' @keywords internal
add_pid_to_tax <- function(tax_results, asv_pid) {
  tax_results$tax <- cbind(tax_results$tax, ASV = rownames(tax_results$tax))
  tax_results$boot <- cbind(tax_results$boot, ASV = asv_pid)
  return(tax_results)
}

#' Combine taxonomic assignments and bootstrap values for each locus into single falsification vector
#'
#' @param tax_results The dataframe containing taxonomic assignments
#' @keywords internal
assignTax_as_char <- function(tax_results) {
  tax_matrix_path <- file.path(directory_path, "Final_tax_matrix.Rdata")
  tax_out <- vapply(1:nrow(tax_results$tax), FUN.VALUE = character(1), function(i) {
    paste(tax_results$tax[i, ],
          tax_results$boot[i, ],
          colnames(tax_results$tax),
          sep = '--', collapse = ';')
  })
  names(tax_out) <- rownames(tax_results$tax)
  save(tax_out, file = tax_matrix_path)
  return(tax_out)
  #check
  stopifnot(all(names(tax_out) %in% colnames(asv_abund_matrix)))
  stopifnot(all(! duplicated(names(tax_out))))
}

#' Format ASV abundance matrix
#'
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#' @param seq_tax_asv An amplified sequence variants matrix with taxonomic information
#' @keywords internal
#Reformat ASV matrix
format_abund_matrix <- function(asv_abund_matrix, seq_tax_asv) {
  formatted_abund_asv <- t(asv_abund_matrix)
  formatted_abund_asv <- cbind(sequence = rownames(formatted_abund_asv),
                               dada2_tax = str_match(seq_tax_asv[rownames(formatted_abund_asv)], pattern = "^(.+)--Species")[,1],
                               dada2_pid = as.numeric(str_match(seq_tax_asv[rownames(formatted_abund_asv)], '--([0-9.]+)--ASV$')[, 2]),
                               formatted_abund_asv)
  formatted_abund_asv <- as_tibble(formatted_abund_asv)
  write_csv(formatted_abund_asv, file = file.path(directory_path, 'final_asv_abundance_matrix.csv'))
  #make results folder
  print(formatted_abund_asv)
  return(formatted_abund_asv)
}

#' Final inventory of read counts after each step from input to removal of chimeras. This function deals with if you have more than one sample. TODO optimize for one sample
#'
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#' @keywords internal
get_read_counts <- function(asv_abund_matrix) {
  filter_results_path<-file.path(directory_path, "filter_results.RData")
  load(filter_results_path) #incorporate into function
  denoised_data_path <- file.path(directory_path, "Denoised_data.Rdata")
  load(denoised_data_path) #incorporate into function
  merged_read_data_path <- file.path(directory_path, "Merged_reads.Rdata")
  load(merged_read_data_path)
  getN <- function(x) sum(dada2::getUniques(x))
  track <- cbind(filter_results, sapply(dada_forward, getN), sapply(dada_reverse, getN), sapply(merged_reads, getN), rowSums(asv_abund_matrix))
  track <- cbind(rownames(track), data.frame(track, row.names=NULL))
  colnames(track) <- c("sample_nameBarcode", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  track$sample_nameBarcode<- gsub(".fastq.gz", "",
                                  gsub("R1_", "", track$sample_nameBarcode, fixed = TRUE))
  track_read_counts_path <- file.path(directory_path, "track_reads.csv")
  write_csv(track, track_read_counts_path)
  print(track)
}
