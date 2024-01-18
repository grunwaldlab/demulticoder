#' Assign taxonomy functions
#' @param analysis_setup A list containing directory paths and data tables, produced by the `prepare_reads` function.
#' @param barcode specify which barcode you have used, 'rps10', 'its', or 'rps10_its'
#' @param asv_abund_matrix specify the ASV abundance matrix for which taxonomic assignments will be given
#' @param retrieve_files Specify TRUE/FALSE whether to copy files from the tempdirectory (which will be deleted) 
#' to directory specified by directory path
#' @param overwrite_existing Logical, indicating whether to remove or overwrite existing files and directories from previous runs. If set to TRUE, specific output files 
#' @inheritParams assign_taxonomyDada2
#' @inheritParams dada2::assignTaxonomy
#' @return Taxonomic assignments of each unique ASV sequence
#' @export assignTax
#' @examples
#' Assign taxonomies to ASVs on a per barcode basis
#' prepare_reads(maxN = 0, data_directory = "~/demulticoder/inst/extdata", output_directory = "~/testing_package", tempdir_id = "run1", overwrite_existing = TRUE)
#' cut_trim(analysis_setup,cutadapt_path="/opt/homebrew/bin/cutadapt", overwrite_existing = TRUE)
#' make_asv_abund_matrix(analysis_setup, overwrite_existing = TRUE)
#' assignTax(analysis_setup,asv_abund_matrix, retrieve_files=TRUE, overwrite_existing=TRUE)

assignTax <- function(analysis_setup, asv_abund_matrix, tryRC = FALSE, verbose = FALSE, multithread = FALSE, retrieve_files = FALSE, db_rps10="oomycetedb.fasta", db_its="fungidb.fasta", db_16s="bacteriadb.fasta", db_other="otherdb.fasta", overwrite_existing=FALSE) {
  dir_paths <- analysis_setup$dir_paths
  data_tables <- analysis_setup$data_tables
  directory_path <- dir_paths$output_directory
  data_path <- dir_paths$data_directory
  directory_path_temp <- dir_paths$temp_directory
  unique_barcodes <- unique(data_tables$cutadapt_data$primer_name)
  
  if (overwrite_existing) {
    
    patterns_to_remove <- c(
      "*_species_count_table.csv",
      "dada2_asv_alignments.txt",
      "final_asv_abundance_matrix.csv",
      "track_reads.csv"
    )
    
    for (pattern in patterns_to_remove) {
      full_pattern <- file.path(directory_path, pattern)
      files_to_remove <- list.files(path = directory_path, pattern = pattern, full.names = TRUE)
      
      if (length(files_to_remove) > 0) {
        file.remove(files_to_remove)
      }
    }
    
    patterns_to_remove_temp <- c(
      "*_reference_db.fa",
      "*_taxmatrix.Rdata",
      "Final_tax_matrix.Rdata"
    )
    
    for (pattern in patterns_to_remove_temp) {
      full_pattern <- file.path(directory_path_temp, pattern)
      files_to_remove <- list.files(path = directory_path_temp, pattern = pattern, full.names = TRUE)
      
      if (length(files_to_remove) > 0) {
        file.remove(files_to_remove)
      }
    }
    
    for (barcode in unique_barcodes) {
      # Load merged reads for the current barcode
      load(file.path(directory_path_temp, paste0("asvabund_matrixDADA2_", barcode, ".Rdata")))
      
      format_database(analysis_setup, barcode, db_its, db_rps10, db_16s, db_other)
      
      # Run taxonomy assignment for the current barcode
      process_single_barcode(
        data_tables=data_tables,
        directory_path_temp = directory_path_temp,
        directory_path = dir_paths$output_directory,
        asv_abund_matrix = asv_abund_matrix,
        locus = barcode
      )
    }
    if (retrieve_files) {
      file.copy(directory_path_temp, dir_paths$output_directory, recursive = TRUE)
    }
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
           directory_path_temp,
           directory_path,
           asv_abund_matrix,
           tryRC = FALSE,
           verbose = FALSE,
           multithread = FALSE,
           locus = barcode)
  {
    abund_asv_single <-
      prep_abund_matrix(data_tables$cutadapt_data, asv_abund_matrix, data_tables, locus)
    refdb = paste0(locus, "_reference_db.fa")
    taxmat = paste0(locus, "_taxmatrix.Rdata")
    tax_results_single_asv <-
      assign_taxonomyDada2(
        abund_asv_single,
        directory_path_temp,
        tryRC = tryRC,
        verbose = verbose,
        multithread = multithread,
        locus=locus
      )
    single_pids_asv <- get_pids(tax_results_single_asv, directory_path_temp, directory_path, refdb, locus)
    tax_results_single_asv_pid <-
      add_pid_to_tax(tax_results_single_asv, single_pids_asv)
    seq_tax_asv <- assignTax_as_char(tax_results_single_asv_pid, directory_path_temp, locus)
    formatted_abund_asv <-
      format_abund_matrix(asv_abund_matrix, seq_tax_asv, directory_path, locus)
    get_read_counts(asv_abund_matrix, directory_path_temp, directory_path, locus)
  }

###
#' Prepare final ASV abundance matrix
#'
#' @param directory_data folder with trimmed and filtered reads for each sample
#' @param asv_abund_matrix The returned final ASV abundance matrix
#' @param locus The barcode selected in the analysis
#' @keywords internal
prep_abund_matrix <-function(cutadapt_data, asv_abund_matrix, data_tables, locus){
  #rownames(asv_abund_matrix) <- sub(rownames(asv_abund_matrix), pattern = ".fastq.gz", replacement = "")
  data_tables$cutadapt_data$file_id_primer <- paste0(data_tables$cutadapt_data$sample_name, "_", data_tables$cutadapt_data$primer_name)
  asv_abund_matrix<- asv_abund_matrix[rownames(asv_abund_matrix) %in% data_tables$cutadapt_data$file_id_primer[data_tables$cutadapt_data$primer_name == locus], ]
  return(asv_abund_matrix)
}


#' Assign taxonomy
#'
#' @inheritParams dada2::assignTaxonomy
#' @param asv_abund_matrix The ASV abundance matrix
#' @param ref_database The reference database used for taxonomic inference steps
#' @keywords internal
assign_taxonomyDada2<-function(asv_abund_matrix, directory_path_temp, minBoot=0, tryRC=FALSE, verbose=FALSE, multithread=TRUE, locus=barcode){
  tax_results<- dada2::assignTaxonomy(asv_abund_matrix,
                                      refFasta = file.path(directory_path_temp, paste0(locus, "_reference_db.fa")),
                                      taxLevels = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Reference"),
                                      minBoot = minBoot,
                                      tryRC = tryRC,
                                      outputBootstraps = TRUE,
                                      multithread = multithread)
  tax_matrix_path <- file.path(directory_path_temp, paste0("TaxMatrix_", locus, ".Rdata"))
  save(tax_results, file = tax_matrix_path)
  return(tax_results)
}

##

#' Align ASV sequences to reference sequences from database to get percent ID. Get percent identities.
#'
#' @param tax_results The dataframe containing taxonomic assignments
#' @keywords internal

get_pids <- function(tax_results, directory_path_temp, directory_path, db, locus) {
  db_seqs <- read_fasta(file.path(directory_path_temp, db))
  ref_seqs <- get_ref_seq(tax_results, db_seqs)
  asv_seqs <- rownames(tax_results$tax)
  
  # Calculate pids for normal alignment
  asv_ref_align_normal <- mapply(pairwiseAlignment, asv_seqs, ref_seqs, MoreArgs = list(type = 'global-local'))
  asv_pids_normal <- sapply(asv_ref_align_normal, function(align) {
    pid(align, type = "PID1")
  })
  
  # Calculate pids for reverse complement alignment
  asv_ref_align_revcomp <- mapply(pairwiseAlignment, asv_seqs, rev_comp(ref_seqs), MoreArgs = list(type = 'global-local'))
  asv_pids_revcomp <- sapply(asv_ref_align_revcomp, function(align) {
    pid(align, type = "PID1")
  })
  
  # Choose the alignment with higher pid
  asv_ref_align <- Map(function(align_normal, align_revcomp, pid_normal, pid_revcomp) {
    if (pid_normal >= pid_revcomp) align_normal else align_revcomp
  }, asv_ref_align_normal, asv_ref_align_revcomp, asv_pids_normal, asv_pids_revcomp)
  
  # Choose the higher pid value
  asv_pids <- mapply(function(pid_normal, pid_revcomp) {
    if (pid_normal >= pid_revcomp) pid_normal else pid_revcomp
  }, asv_pids_normal, asv_pids_revcomp)
  
  # Write the alignment results to a text file
  alignment_text <- paste0(
    'asv: ', asv_seqs, '\n',
    'ref: ', unlist(ref_seqs), '\n',
    'pid: ', asv_pids, '\n',
    unlist(sapply(asv_ref_align, function(align) align@pattern)), '\n',
    unlist(sapply(asv_ref_align, function(align) align@subject)), '\n'
  )
  write_lines(paste(collapse = '==================================================\n', alignment_text),
              file = file.path(directory_path, paste0("dada2_asv_alignments_", locus, ".txt")))
  
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
assignTax_as_char <- function(tax_results, directory_path_temp, locus) {
  tax_matrix_path <- file.path(directory_path_temp, paste0("Final_tax_matrix_", locus, ".Rdata"))
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
format_abund_matrix <- function(asv_abund_matrix, seq_tax_asv, directory_path, locus) {
  formatted_abund_asv <- t(asv_abund_matrix)
  formatted_abund_asv <- cbind(sequence = rownames(formatted_abund_asv),
                               dada2_tax = str_match(seq_tax_asv[rownames(formatted_abund_asv)], pattern = "^(.+)--Species")[,1],
                               dada2_pid = as.numeric(str_match(seq_tax_asv[rownames(formatted_abund_asv)], '--([0-9.]+)--ASV$')[, 2]),
                               formatted_abund_asv)
  formatted_abund_asv <- as_tibble(formatted_abund_asv)
  
  primer_seqs <- apply(analysis_setup$data_tables$primer_data[, 2:ncol(analysis_setup$data_tables$primer_data)], 2, paste, collapse = "|")
  
  # Function to create a regular expression pattern for a primer with ambiguity codes
  # Function to create a regular expression pattern for a primer with ambiguity codes
  create_primer_pattern <- function(primer_seq) {
    # Replace ambiguity codes with corresponding regular expression patterns
    ambig_code_mapping <- c("R" = "[AG]", "Y" = "[CT]", "S" = "[GC]", "W" = "[AT]", "K" = "[GT]", "M" = "[AC]", "B" = "[CGT]", "D" = "[AGT]", "H" = "[ACT]", "V" = "[ACG]", "N" = "[ACGT]")
    
    pattern <- ""
    for (char in strsplit(primer_seq, '')[[1]]) {
      if (char %in% names(ambig_code_mapping)) {
        pattern <- paste0(pattern, ambig_code_mapping[char])
      } else {
        pattern <- paste0(pattern, char)
      }
    }
    
    return(pattern)
  }
  
  # Create regular expression patterns for each primer
  primer_patterns <- sapply(primer_seqs, create_primer_pattern)
  
  # Filter out sequences based on primers
  keep_rows <- sapply(formatted_abund_asv$sequence, function(asv_seq) {
    # Check if any of the primer sequences are present in the ASV sequence
    match_primer <- sapply(primer_patterns, function(pattern) any(grepl(pattern, asv_seq, ignore.case = TRUE)))
    
    if (any(match_primer)) {
      cat("ASV Sequence:", asv_seq, "\n")
      cat("Matching Primer(s):", primer_seqs[match_primer], "\n")
      cat("Sequence removed\n-----\n")
    }
    
    # Keep the row if none of the primer sequences match
    !any(match_primer)
  })
  
  # Update the abundance matrix with filtered rows
  filtered_abund_matrix <- formatted_abund_asv[keep_rows, ]
  
  # Return the filtered abundance matrix
  
  write_csv(filtered_abund_matrix, file = file.path(directory_path, paste0('final_asv_abundance_matrix_', locus, '.csv')))
  return(filtered_abund_matrix)
}


#' Final inventory of read counts after each step from input to removal of chimeras. This function deals with if you have more than one sample. TODO optimize for one sample
#'
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#' @keywords internal
get_read_counts <- function(asv_abund_matrix, directory_path_temp, directory_path, locus) {
  filter_results_path<-file.path(directory_path_temp, paste0("Filter_results_", locus, ".RData"))
  load(filter_results_path) #incorporate into function
  denoised_data_path <- file.path(directory_path_temp, paste0("Denoised_data_", locus, ".RData"))
  load(denoised_data_path) #incorporate into function
  merged_read_data_path <- file.path(directory_path_temp, paste0("Merged_reads_", locus, ".Rdata"))
  load(merged_read_data_path)
  getN <- function(x) sum(dada2::getUniques(x))
  track <- cbind(filter_results, sapply(dada_forward, getN), sapply(dada_reverse, getN), sapply(merged_reads, getN), rowSums(asv_abund_matrix))
  track <- cbind(rownames(track), data.frame(track, row.names=NULL))
  colnames(track) <- c("samplename_barcode", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  track$samplename_barcode<- gsub(".fastq.gz", "",
                                  gsub("R1_", "", track$samplename_barcode, fixed = TRUE))
  track_read_counts_path <- file.path(directory_path, paste0("track_reads_", locus, ".csv"))
  write_csv(track, track_read_counts_path)
  print(track)
}