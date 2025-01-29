#' Prepare final ASV abundance matrix
#'
#' @param directory_data folder with trimmed and filtered reads for each sample
#' @param asv_abund_matrix The returned final ASV abundance matrix
#' @param locus The barcode selected in the analysis
#'
#' @keywords internal
prep_abund_matrix <-function(cutadapt_data, asv_abund_matrix, data_tables, locus){
  data_tables$cutadapt_data$file_id_primer <- paste0(data_tables$cutadapt_data$sample_name, "_", data_tables$cutadapt_data$primer_name)
  asv_abund_matrix<- asv_abund_matrix[rownames(asv_abund_matrix) %in% data_tables$cutadapt_data$file_id_primer[data_tables$cutadapt_data$primer_name == locus], ]
  return(asv_abund_matrix)
}
#' Assign taxonomy
#'
#' @inheritParams dada2::assignTaxonomy
#' @inheritParams RcppParallel::setThreadOptions
#' @param asv_abund_matrix The ASV abundance matrix
#' @param ref_database The reference database used for taxonomic inference steps
#' 
#' @keywords internal
assign_taxonomyDada2<-function(asv_abund_matrix, temp_directory_path, minBoot=0, tryRC=FALSE, verbose=FALSE, multithread=TRUE, locus=barcode){
  set.seed(1) #add parameter
  tax_results<- dada2::assignTaxonomy(asv_abund_matrix,
                                      refFasta = file.path(temp_directory_path, paste0(locus, "_reference_db.fa")),
                                      taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Reference"),
                                      minBoot = minBoot,
                                      tryRC = tryRC,
                                      outputBootstraps = TRUE,
                                      multithread = multithread)
  tax_matrix_path <- file.path(temp_directory_path, paste0("TaxMatrix_", locus, ".RData"))
  save(tax_results, file = tax_matrix_path)
  return(tax_results)
}

#' Align ASV sequences to reference sequences from database to get percent ID.
#' Get percent identities.
#'
#' @param tax_results The data frame containing taxonomic assignments
#'
#' @keywords internal
get_pids <- function(tax_results, temp_directory_path, output_directory_path, db, locus) {
  db_seqs <- metacoder::read_fasta(file.path(temp_directory_path, db))
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
  # Choose the higher pid value
  asv_pids <- ifelse(asv_pids_normal >= asv_pids_revcomp, asv_pids_normal, asv_pids_revcomp)
  # Format the alignment results
  alignment_text <- paste0(map_chr(seq_along(asv_ref_align_normal), function(i) {
    paste0('ASV Sequence: ', asv_seqs[[i]], '\n',
           'Reference Sequence: ', unlist(ref_seqs[[i]]),'\n',
           'Percent Identity: ', asv_pids[[i]],'\n',
           'Aligned ASV Sequence: ', asv_ref_align_normal[[i]]@pattern,'\n',
           'Aligned Reference Sequence: ',asv_ref_align_normal[[i]]@subject,'\n')
  }), 
  collapse = '==================================================\n')
  # Write the alignment results to a text file
  readr::write_lines(alignment_text, file = file.path(output_directory_path, paste0("dada2_asv_alignments_", locus, ".txt")))
  return(asv_pids)
}

#' Align ASV sequences to reference sequences from database to get percent ID.
#' Start by retrieving reference sequences.
#'
#' @param tax_results The dataframe containing taxonomic assignments
#' @param db The reference database
#'
#' @keywords internal
get_ref_seq <- function(tax_results, db) {
  ref_i <- as.integer(stringr_str_match(tax_results$tax[, 'Reference'], '^.+_([0-9]+)$')[ ,2])
  db[ref_i]
}

#' Add PID and bootstrap values to tax result.
#'
#' @param tax_results The dataframe containing taxonomic assignments
#' @param asv_pid Percent identity information for each ASV relative to
#'   reference database sequence
#' @keywords internal
add_pid_to_tax <- function(tax_results, asv_pid) {
  tax_results$tax <- cbind(tax_results$tax, ASV = rownames(tax_results$tax))
  tax_results$boot <- cbind(tax_results$boot, ASV = asv_pid)
  return(tax_results)
}

#' Combine taxonomic assignments and bootstrap values for each locus into single
#' falsification vector
#'
#' @param tax_results The dataframe containing taxonomic assignments
#'
#' @keywords internal
assignTax_as_char <- function(tax_results, temp_directory_path, locus) {
  tax_matrix_path <- file.path(temp_directory_path, paste0("Final_tax_matrix_", locus, ".RData"))
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
#' @param data_tables The data tables containing the paths to read files, metadata, primer sequences
#' @param asv_abund_matrix An abundance matrix containing amplified sequence
#'   variants
#' @param seq_tax_asv An amplified sequence variants matrix with taxonomic
#'   information
#'
#' @keywords internal
format_abund_matrix <- function(data_tables, asv_abund_matrix, seq_tax_asv, output_directory_path, locus) {
  formatted_abund_asv <- t(asv_abund_matrix)
  #asv_id_column <- paste("asv_", seq_along(rownames(formatted_abund_asv)), sep = "")
  formatted_abund_asv <- cbind(
    #asv_id = asv_id_column,
    sequence = rownames(formatted_abund_asv),
    dada2_tax = stringr::str_match(seq_tax_asv[rownames(formatted_abund_asv)], pattern = "^(.+)--Species")[,1],
    #dada2_pid = as.numeric(stringr_str_match(seq_tax_asv[rownames(formatted_abund_asv)], '--([0-9.]+)--ASV$')[, 2]),
    formatted_abund_asv)
  formatted_abund_asv <- tibble::as_tibble(formatted_abund_asv)
  
  primer_seqs <- apply(data_tables$primer_data[, 2:ncol(data_tables$primer_data)], 2, paste, collapse = "|")
  create_primer_pattern <- function(primer_seq) {
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
  asv_id_column <- paste("asv_", seq_along(rownames(filtered_abund_matrix)), sep = "")
  filtered_abund_matrix <- cbind(
    asv_id = asv_id_column,
    filtered_abund_matrix)

  # Return the filtered abundance matrix
  
  readr::write_csv(filtered_abund_matrix, file = file.path(output_directory_path, paste0('final_asv_abundance_matrix_', locus, '.csv')))
  return(filtered_abund_matrix)
}

#' Final inventory of read counts after each step from input to removal of chimeras. This function deals with if you have more than one sample. TODO optimize for one sample
#'
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#' 
#' @keywords internal
get_read_counts <- function(asv_abund_matrix, temp_directory_path, output_directory_path, locus) {
  filter_results_path <- file.path(temp_directory_path, paste0("Filter_results_", locus, ".RData"))
  load(filter_results_path)
  denoised_data_path <- file.path(temp_directory_path, paste0("Denoised_data_", locus, ".RData"))
  load(denoised_data_path)
  merged_read_data_path <- file.path(temp_directory_path, paste0("Merged_reads_", locus, ".RData"))
  load(merged_read_data_path)
  getN <- function(x) sum(dada2::getUniques(x))
  track <- cbind(filter_results, sapply(dada_forward, getN), sapply(dada_reverse, getN), sapply(merged_reads, getN), rowSums(asv_abund_matrix))
  track <- cbind(rownames(track), data.frame(track, row.names=NULL))
  colnames(track) <- c("samplename_barcode", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  track$samplename_barcode <- gsub(".fastq.gz", "",
                                   gsub("R1_", "", track$samplename_barcode, fixed = TRUE))
  
  output_format <- knitr::opts_knit$get("rmarkdown.pandoc.to")
  print(track)
  track_read_counts_path <- file.path(output_directory_path, paste0("track_reads_", locus, ".csv"))
  readr::write_csv(track, track_read_counts_path)
}

#' Process the information from an ASV abundance matrix to run DADA2 for single
#' barcode
#'
#' @param data_tables The data tables containing the paths to read files, metadata, primer sequences
#' @param asv_abund_matrix An abundance matrix containing amplified sequence
#'   variants
#' @inheritParams assign_taxonomyDada2
#' @inheritParams dada2::assignTaxonomy
#' @inheritParams RcppParallel::setThreadOptions
#'
#' @keywords internal
process_single_barcode <-
  function(data_tables,
           temp_directory_path,
           output_directory_path,
           asv_abund_matrix,
           tryRC = FALSE,
           verbose = FALSE,
           multithread = FALSE,
           locus = barcode)
  {
    abund_asv_single <-
      prep_abund_matrix(data_tables$cutadapt_data, asv_abund_matrix, data_tables, locus)
    refdb = paste0(locus, "_reference_db.fa")
    taxmat = paste0(locus, "_taxmatrix.RData")
    tax_results_single_asv <-
      assign_taxonomyDada2(
        abund_asv_single,
        temp_directory_path,
        tryRC = tryRC,
        verbose = verbose,
        multithread = multithread,
        locus=locus
      )
    #single_pids_asv <- get_pids(tax_results_single_asv, temp_directory_path, output_directory_path, refdb, locus)
    #tax_results_single_asv_pid <-
    #add_pid_to_tax(tax_results_single_asv, single_pids_asv)
    seq_tax_asv <- assignTax_as_char(tax_results_single_asv, temp_directory_path, locus)
    formatted_abund_asv <-
      format_abund_matrix(data_tables, asv_abund_matrix, seq_tax_asv, output_directory_path, locus)
    get_read_counts(asv_abund_matrix, temp_directory_path, output_directory_path, locus)
  }


#' Assign taxonomy functions
#' @importFrom utils modifyList read.table stack
#' @param analysis_setup An object containing directory paths and data tables,
#'   produced by the `prepare_reads` function
#' @param asv_abund_matrix ASV abundance matrix.
#' @param tryRC Whether to try reverse complementing sequences during taxonomic
#'   assignment
#' @param db_rps10 The reference database for the rps10 locus
#' @param db_its The reference database for the ITS locus
#' @param db_16S The reference database for the 16S locus
#' @param db_other1 The reference database for different locus 1 (assumes format is like SILVA DB entries)
#' @param db_other2 The reference database for a different locus 2 (assumes format is like SILVA DB entries)
#' @param verbose Logical, indicating whether to display verbose output
#' @param multithread Logical, indicating whether to use multithreading
#' @param retrieve_files Specify TRUE/FALSE whether to copy files from the temp
#'   directory to the output directory
#' @param overwrite_existing Logical, indicating whether to remove or overwrite
#'   existing files and directories from previous runs. Default is `FALSE`.
#'
#' @return Taxonomic assignments of each unique ASV sequence
#'
#' @export assign_tax
#'
#' @examples
#' # Assign taxonomies to ASVs on a per barcode basis
#' analysis_setup <- prepare_reads(
#'   data_directory = system.file("extdata", package = "demulticoder"),
#'   output_directory = tempdir(),
#'   tempdir_path = tempdir(),
#'   tempdir_id = "demulticoder_run_temp",
#'   overwrite_existing = TRUE
#' )
#' cut_trim(
#' analysis_setup,
#' cutadapt_path="/usr/bin/cutadapt",
#' overwrite_existing = TRUE
#' )
#' make_asv_abund_matrix(
#' analysis_setup, 
#' overwrite_existing = TRUE
#' )
#' assign_tax(
#' analysis_setup,
#' asv_abund_matrix, 
#' retrieve_files=FALSE, 
#' overwrite_existing = TRUE
#' )
assign_tax <- function(analysis_setup, asv_abund_matrix, tryRC = FALSE, verbose = FALSE, multithread = FALSE, retrieve_files = FALSE, overwrite_existing = FALSE, db_rps10 = "oomycetedb.fasta", db_its = "fungidb.fasta", db_16S = "bacteriadb.fasta", db_other1 = "otherdb1.fasta", db_other2 = "otherdb2.fasta") {
  data_tables <- analysis_setup$data_tables
  data_path <- analysis_setup$directory_paths$data_directory
  output_directory_path <- analysis_setup$directory_paths$output_directory
  temp_directory_path <- analysis_setup$directory_paths$temp_directory
  unique_barcodes <- unique(data_tables$cutadapt_data$primer_name)
  
  files_to_check <- c("*_reference_db.fa", "TaxMatrix_*", "Final_tax_matrix_*")
  existing_files <- list.files(temp_directory_path, pattern = files_to_check, full.names = TRUE)
  
  if (!overwrite_existing && length(existing_files) > 0) {
    message("Existing files found. Specify overwrite=TRUE, to rerun analysis")
    return(invisible())
    
  } else {
    patterns_to_remove <- c(
      "*_species_count_table.csv",
      "dada2_asv_alignments_*",
      "final_asv_abundance_matrix_*",
      "track_reads_*"
    )
    
    for (pattern in patterns_to_remove) {
      full_pattern <- file.path(output_directory_path, pattern)
      files_to_remove <- list.files(path = output_directory_path, pattern = pattern, full.names = TRUE)
      
      if (length(files_to_remove) > 0) {
        file.remove(files_to_remove)
      }
    }
    
    patterns_to_remove_temp <- c(
      "*_reference_db.fa",
      "TaxMatrix_*",
      "Final_tax_matrix_*.RData"
    )
    
    for (pattern in patterns_to_remove_temp) {
      full_pattern <- file.path(temp_directory_path, pattern)
      files_to_remove <- list.files(path = temp_directory_path, pattern = pattern, full.names = TRUE)
      
      if (length(files_to_remove) > 0) {
        file.remove(files_to_remove)
      }
    }
    
    if (length(existing_files) == 0 && !overwrite_existing) {
      warning("No existing files found. The analysis will be run.")
    }
    
    for (barcode in unique_barcodes) {
      # Load merged reads for the current barcode
      load(file.path(temp_directory_path, paste0("asvabund_matrixDADA2_", barcode, ".RData")))
      
      format_database(data_tables, data_path, output_directory_path, temp_directory_path, barcode, db_its, db_rps10, db_16S, db_other1, db_other2)
      
      # Run taxonomy assignment for the current barcode
      process_single_barcode(
        data_tables = data_tables,
        temp_directory_path = temp_directory_path,
        output_directory_path = output_directory_path,
        asv_abund_matrix = asv_abund_matrix,
        locus = barcode
      )
    }
    
    if (retrieve_files) {
      file.copy(temp_directory_path, output_directory_path,  recursive = TRUE)
    }
  }
  return(invisible())
}