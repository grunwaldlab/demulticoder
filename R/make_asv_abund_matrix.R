#' Make an amplified sequence variant (ASV) abundance matrix This function
#' generates an ASV abundance matrix using raw reads processed during previous
#' steps, including read preparation, removing primers, and using DADA2 core
#' denoising alogrithm to infer ASVs.
#'
#' @param analysis_setup analysis_setup An object containing directory paths and
#'   data tables, produced by the `prepare_reads` function
#' @param overwrite_existing Logical, indicating whether to overwrite existing
#'   results
#' @param ... Additional parameters to pass to dada function
#' @details The function processes data for each unique barcode separately,
#'   inferring ASVs, merging reads, and creating an ASV abundance matrix 
#' @return The ASV abundance matrix (`asv_abund_matrix`)
#' @export make_asv_abund_matrix
#'
#' @examples
#' # The primary wrapper function for DADA2 ASV inference steps
#' prepare_reads(data_directory = "extdata", 
#' output_directory = "test_data", 
#' tempdir_id = "demulticoder_run", 
#' overwrite_existing = FALSE)
#' cut_trim(analysis_setup,
#' cutadapt_path="/opt/homebrew/bin/cutadapt", 
#' overwrite_existing = FALSE)
#' make_asv_abund_matrix(analysis_setup, 
#' overwrite_existing = FALSE)

make_asv_abund_matrix <- function(analysis_setup, overwrite_existing = FALSE) {
  
  data_tables <- analysis_setup$data_tables
  output_directory_path <- analysis_setup$directory_paths$output_directory
  temp_directory_path <- analysis_setup$directory_paths$temp_directory
  
  asv_abund_matrix_list <- list()  
  
  default_params <- list(
    multithread = FALSE,
    nbases = 1e+08,
    randomize = FALSE,
    MAX_CONSIST = 10,
    OMEGA_C = 0,
    qualityType = "Auto",
    nominalQ = FALSE,
    obs = TRUE,
    err_out = TRUE,
    err_in = FALSE,
    pool = FALSE,
    selfConsist = FALSE,
    verbose = FALSE,
    minOverlap = 12,
    maxMismatch = 0,
    method = "consensus",
    min_asv_length = 50)
  
  files_to_check <- c("asvabund_matrixDADA2_*")
  existing_files <- list.files(temp_directory_path, pattern = files_to_check, full.names = TRUE)
  
  if (!overwrite_existing && length(existing_files) > 0) {
    message("Existing files found. The following ASV abundance matrices are saved in the tempdir path:")
    
    asv_abund_matrix_list <- lapply(existing_files, function(file) {
      message(file)
    })
    return(asv_abund_matrix_list)
  } else {
    patterns_to_remove <- c(
      "asv_seqlength_plot_*",
      "error_plot_*",
      "read_merging_info_*"
    )
    
    for (pattern in patterns_to_remove) {
      full_pattern <- file.path(output_directory_path, pattern)
      files_to_remove <- list.files(path = output_directory_path, pattern = pattern, full.names = TRUE)
      
      if (length(files_to_remove) > 0) {
        file.remove(files_to_remove)
      }
    }
    
    patterns_to_remove_temp <- "asvabund_matrixDADA2_*"
    
    for (pattern in patterns_to_remove_temp) {
      full_pattern <- file.path(temp_directory_path, pattern)
      files_to_remove <- list.files(path = temp_directory_path, pattern = pattern, full.names = TRUE)
      
      if (length(files_to_remove) > 0) {
        file.remove(files_to_remove)
      }
    }
    
    if (length(existing_files) == 0 && !overwrite_existing) {
      warning("No existing files found. The 'make_asv_abund_matrix' function will run.")
      # Continue with the analysis
    }
    
    unique_barcodes <- unique(data_tables$cutadapt_data$primer_name)
    
    for (barcode in unique_barcodes) {
      # Get barcode-specific parameters
      params <- dplyr::filter(data_tables$parameters, primer_name == barcode)
      if (nrow(params) > 0) {
        params <- as.list(params)
        
        barcode_params <- modifyList(default_params, params)
        
        if (overwrite_existing || !file.exists(file.path(temp_directory_path, paste0("asvabund_matrixDADA2_", barcode, ".RData")))) {
          infer_asv_command(
            output_directory_path = output_directory_path,
            temp_directory_path = temp_directory_path,
            data_tables = data_tables,
            barcode_params = barcode_params,
            barcode = barcode
          )
          
          # Merge reads with barcode-specific parameters
          merged_reads <- merge_reads_command(
            output_directory_path = output_directory_path,
            temp_directory_path = temp_directory_path,
            barcode_params = barcode_params,
            barcode = barcode
          )
          
          countOverlap(merged_reads, data_tables, output_directory_path, barcode)
          raw_seqtab <- createASVSequenceTable(merged_reads, orderBy = "abundance")
          asv_abund_matrix <- make_abund_matrix(raw_seqtab, temp_directory_path = temp_directory_path, barcode_params, barcode)
          make_seqhist(asv_abund_matrix, output_directory_path)
          asv_abund_matrix_list[[barcode]] <- file.path(temp_directory_path, paste0("asvabund_matrixDADA2_", barcode, ".RData"))
        }
      }
    }
    
    return(asv_abund_matrix_list)
  }
}

#' Retrieve the paths of the filtered and trimmed Fastq files
#'
#' @param my_direction Whether primer is in forward or reverse direction
#' @param my_primer_pair_id The the specific barcode id
#' @param cutadapt_data directory_data folder with trimmed and filtered reads
#'   for each sample
#' @keywords internal
get_fastq_paths <- function(analysis_setup, my_direction, my_primer_pair_id) {
  data_tables <- analysis_setup$data_tables
  filtered_paths <- character()
  for (i in seq_along(data_tables$cutadapt_data$direction)) {
    if (data_tables$cutadapt_data$direction[i] == my_direction &&
        data_tables$cutadapt_data$primer_name[i] == my_primer_pair_id &&
        file.exists(data_tables$cutadapt_data$filtered_path[i])) {
      filtered_paths <- c(filtered_paths, data_tables$cutadapt_data$filtered_path[i])
    }
  }
  
  filtered_paths
}
#' Core DADA2 function to learn errors and infer ASVs
#' @param output_directory_path The path to the directory containing the fastq,
#' metadata, and primerinfo_params files
#' @param my_primer_pair_id The the specific barcode id 
#' @param my_direction Location of read files and metadata file
#' @return asv_data
#' @keywords internal
infer_asvs <- function(my_direction, my_primer_pair_id, barcode_params, output_directory_path) {
  set.seed(1)
  fastq_paths <- get_fastq_paths(analysis_setup, my_direction, my_primer_pair_id)
  error_plot_filename <- paste0("error_plot_", my_primer_pair_id, ".pdf")
  error_profile <- dada2::learnErrors(
    fastq_paths,
    multithread = barcode_params$multithread,
    nbases = barcode_params$nbases,
    randomize = barcode_params$randomize,
    MAX_CONSIST = barcode_params$MAX_CONSIST,
    OMEGA_C = barcode_params$OMEGA_C,
    qualityType = barcode_params$qualityType,
    verbose = barcode_params$verbose
  )
  
  cat(
    paste0(
      'Error rate plot for the ',
      my_direction,
      ' read of primer pair ',
      my_primer_pair_id,
      ' \n'
    )
  )
  
  plot_errors <- dada2::plotErrors(
    error_profile,
    nominalQ = barcode_params$nominalQ,
    obs = barcode_params$obs,
    err_out = barcode_params$err_out,
    err_in = barcode_params$err_in
  )
  
  ggplot2::ggsave(
    plot_errors,
    filename = error_plot_filename,
    path = output_directory_path,
    width = 8,
    height = 8
  )
  
  asv_data <- dada2::dada(
    fastq_paths,
    err = error_profile,
    multithread = barcode_params$multithread
  )
  
  return(asv_data)
}

#' Function to infer ASVs, for multiple loci
#'
#' @param output_directory_path The path to the directory where resulting files are output
#' @param denoised_data_path Path to saved intermediate denoised data
#' @keywords internal
infer_asv_command <- function(output_directory_path, temp_directory_path, data_tables, barcode_params, barcode) {
  
  multithread <- barcode_params$multithread
  nbases <- barcode_params$nbases
  randomize <- barcode_params$randomize
  MAX_CONSIST <- barcode_params$MAX_CONSIST
  OMEGA_C <- barcode_params$OMEGA_C
  qualityType <- barcode_params$qualityType
  nominalQ <- barcode_params$nominalQ
  obs <- barcode_params$obs
  err_out <- barcode_params$err_out
  err_in <- barcode_params$err_in
  pool <- barcode_params$pool
  selfConsist <- barcode_params$selfConsist
  verbose <- barcode_params$verbose
  
  denoised_data_path <- file.path(temp_directory_path, paste0("Denoised_data_", barcode, ".RData"))
  
  run_dada <- function(direction, data_tables, barcode_params, barcode) {
    dada_output <- lapply(barcode, function(primer_name) {
      infer_asvs(
        direction,
        primer_name,
        barcode_params, 
        output_directory_path
      )
    })
    unlist(dada_output, recursive = FALSE)
  }
  
  dada_forward <- run_dada("Forward", data_tables, barcode_params, barcode)
  dada_reverse <- run_dada("Reverse", data_tables, barcode_params, barcode)
  save(dada_forward, dada_reverse, file = denoised_data_path)
}


#' Merge forward and reverse reads
#'
#' @inheritParams dada2::mergePairs
#' @param output_directory_path The path to the directory where resulting files
#'   are output
#' @param merged_read_data_path Path to R data file containing merged read data
#' @return merged_reads Intermediate merged read R data file
#' @keywords internal
merge_reads_command <- function(output_directory_path, temp_directory_path, barcode_params, barcode) {
  denoised_data_path <- file.path(temp_directory_path, paste0("Denoised_data_", barcode, ".RData"))
  load(denoised_data_path)
  
  merged_read_data_path <- file.path(temp_directory_path, paste0("Merged_reads_", barcode, ".RData"))
  
  merged_reads <- dada2::mergePairs(
    dadaF = dada_forward,
    derepF = file.path(temp_directory_path, 'filtered_sequences', names(dada_forward)),
    dadaR = dada_reverse,
    derepR = file.path(temp_directory_path, 'filtered_sequences', names(dada_reverse)),
    minOverlap = barcode_params$minOverlap,
    maxMismatch = barcode_params$maxMismatch,
    returnRejects = FALSE,
    justConcatenate = FALSE,
    verbose = barcode_params$verbose
  )
  
  names(merged_reads) <- gsub(".fastq.gz", "", gsub("R1_", "", names(merged_reads), fixed = TRUE))
  save(merged_reads, file = merged_read_data_path)
  return(merged_reads)
}

#' Count overlap to see how well the reads were merged
#'
#' @param merged_reads Intermediate merged read R data file
#' @return A plot describing how well reads merged and information on overlap
#'   between reads Count overlap to see how well the reads were merged
#'
#' @param merged_reads Intermediate merged read R data file
#' @param output_directory_path The path to the directory where resulting files
#'   are output
#' @return A plot describing how well reads merged and information on overlap
#'   between reads
countOverlap <- function(merged_reads, data_tables, output_directory_path, barcode) {
  non_empty_merged_reads <- merged_reads[sapply(merged_reads, nrow) > 0]
  merge_data <- do.call(rbind, non_empty_merged_reads)
  merge_data$samplename_barcode <- rep(names(non_empty_merged_reads), sapply(non_empty_merged_reads, nrow))
  merge_data2 <- merge_data[merge_data$samplename_barcode %in% data_tables$cutadapt_data$samplename_barcode, ]
  merge_data2 <- merge(merge_data2, data_tables$cutadapt_data, by = "samplename_barcode")
  merge_data2$overlap <- merge_data2$nmatch + merge_data2$nmismatch
  merge_data2$mismatch <- merge_data2$nmismatch + merge_data2$nindel
  merge_data2$identity <- (merge_data2$overlap - merge_data2$mismatch) / merge_data2$overlap
  
  # Make dataframe long
  merge_plot <- data.frame(
    barcode = merge_data2$primer_name,
    Mismatches_and_Indels = merge_data2$mismatch,
    Merged = merge_data2$accept,
    Overlap_Length = merge_data2$overlap
  )
  
  long_df <- stack(merge_plot[,c("Mismatches_and_Indels", "Overlap_Length")])
  merge_plot <- cbind(merge_plot[,"barcode"], merge_plot[,"Merged"], long_df)
  names(merge_plot) <- c("barcode", "Merged", "value", "stat")
  
  merge_plot_output <- ggplot2::ggplot(merge_plot, ggplot2::aes(x = value, fill = Merged)) +
    ggplot2::facet_grid(barcode ~ stat, scales = 'free') +
    ggplot2::geom_histogram(bins = 50) +
    ggplot2::scale_fill_viridis_d(begin = 0.8, end = 0.2) +
    ggplot2::labs(x = '', y = 'ASV count', fill = 'Merged') +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
  
  plot_filename <- paste0("read_merging_info_", barcode, ".pdf")
  
  ggplot2::ggsave(
    merge_plot_output,
    filename = plot_filename,
    path = output_directory_path,
    width = 8,
    height = 8
  )
  
  print(merge_plot_output)
}

#' Make ASV sequence matrix
#'
#' @inheritParams dada2::makeSequenceTable
#' @param merged_reads Intermediate merged read R data file
#' @return raw_seqtab
#' @keywords internal
createASVSequenceTable <- function(merged_reads, orderBy = "abundance") {
  raw_seqtab <- dada2::makeSequenceTable(merged_reads, orderBy = orderBy)
  return(raw_seqtab)
}

#' Quality filtering to remove chimeras and short sequences
#'
#' @inheritParams dada2::removeBimeraDenovo
#' @param raw_seqtab R data file with raw sequence data prior to removal of
#'   chimeras
#' @return asv_abund_matrix The returned final ASV abundance matrix
#' @keywords internal
make_abund_matrix <- function(raw_seqtab,
                              temp_directory_path,
                              barcode_params=barcode_params,
                              barcode) {
  seqtab.nochim <- dada2::removeBimeraDenovo(raw_seqtab, method=barcode_params$method, verbose=barcode_params$verbose)
  asv_abund_matrix <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) >= barcode_params$min_asv_length]
  
  asvabund_matrix_path <- file.path(temp_directory_path, paste0("asvabund_matrixDADA2_", barcode, ".RData"))
  save(asv_abund_matrix, file = asvabund_matrix_path)
  return(asv_abund_matrix)
}

#' Plots a histogram of read length counts of all sequences within the ASV
#' matrix
#'
#' @param asv_abund_matrix The returned final ASV abundance matrix
#' @keywords internal
#' @return histogram with read length counts of all sequences within ASV matrix
#' @keywords internal
make_seqhist <- function(asv_abund_matrix, output_directory_path) {
  getPrimerLengths <- function(asv_abund_matrix) {
    barcodes <- unique(sub(".*_(.*)", "\\1", rownames(asv_abund_matrix)))
    primer_lengths <- list()
    
    for (barcode in barcodes) {
      indices <- grepl(paste0("_", barcode), rownames(asv_abund_matrix))
      primer_seqs <- colnames(asv_abund_matrix)[apply(asv_abund_matrix[indices, , drop = FALSE] > 0, 2, any)]
      seq_lengths <- nchar(dada2::getSequences(asv_abund_matrix[indices, primer_seqs, drop = FALSE]))
      primer_lengths[[barcode]] <- seq_lengths
    }
    return(primer_lengths)
  }
  
  primer_lengths <- getPrimerLengths(asv_abund_matrix)
  
  for (barcode in names(primer_lengths)) {
    data <- data.frame(Length = unlist(primer_lengths[[barcode]]))
    
    hist_plot <- ggplot2::ggplot(data, ggplot2::aes(x = Length)) +
      ggplot2::geom_histogram(binwidth = 10, fill = "blue", color = "black", alpha = 0.7, ) +
      ggplot2::labs(x = 'Length of sequence (bp)', y = 'Counts', title = paste("ASV lengths for", barcode, "locus")) +
      ggplot2::theme_minimal()+
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    
    ggplot2::ggsave(
      hist_plot,
      filename = paste("asv_seqlength_plot_", barcode, ".pdf", sep = ""),
      path = output_directory_path,
      width = 8,
      height = 8
    )
    print(hist_plot)
  }
}