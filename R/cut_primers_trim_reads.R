utils::globalVariables(c("f_compt", "f_rc", "f_rev", "forward", "primer_type", "r_compt", "r_rc", "r_rev", "reverse", "filter_results"))

#' Core function for running cutadapt
#'
#' @param cutadapt_path A path to the cutadapt program.
#' @param cutadapt_data FASTQ read files trimmed of primers
#' @param minCutadaptlength Read lengths that are lower than this threshold will
#'   be discarded. Default is 0.
#'
#' @return Trimmed read.
#'
#' @keywords internal
run_cutadapt <- function(cutadapt_path,
                         cutadapt_data_barcode,
                         barcode_params,
                         minCutadaptlength) {
  
  cutadapt <- cutadapt_path
  tryCatch(
    {
      version_output <- system2(cutadapt, args = "--version", stdout = TRUE, stderr = TRUE)
    },
    warning = function(w) {
      stop("cutadapt cannot be found on PATH. Check if it is installed?")
    }
  )
  
  cat("Running Cutadapt", version_output, "for", cutadapt_data_barcode$primer_name[1], "sequence data", "\n")
  
  # Simplify
  cutadapt <- path.expand(cutadapt_path)
  R1_flags = unique(paste("-g", cutadapt_data_barcode$forward, "-a", cutadapt_data_barcode$r_rc))
  R2_flags = unique(paste("-G", cutadapt_data_barcode$reverse, "-A", cutadapt_data_barcode$f_rc))
  fwd_trim = cutadapt_data_barcode[cutadapt_data_barcode$direction == "Forward",][["trimmed_path"]]
  rev_trim = cutadapt_data_barcode[cutadapt_data_barcode$direction == "Reverse",][["trimmed_path"]]
  fwd_untrim = cutadapt_data_barcode[cutadapt_data_barcode$direction == "Forward",][["untrimmed_path"]]
  rev_untrim = cutadapt_data_barcode[cutadapt_data_barcode$direction == "Reverse",][["untrimmed_path"]]
  fwd_prefilt = cutadapt_data_barcode[cutadapt_data_barcode$direction == "Forward",][["prefiltered_path"]]
  rev_prefilt = cutadapt_data_barcode[cutadapt_data_barcode$direction == "Reverse",][["prefiltered_path"]]
  
  # Construct command arguments
  command_args = paste(
    R1_flags,
    R2_flags,
    "-n",
    2,
    "-o",
    fwd_trim,
    "-p",
    rev_trim,
    "--minimum-length",
    minCutadaptlength,
    "--untrimmed-output",
    fwd_untrim,
    "--untrimmed-paired-output",
    rev_untrim,
    "--quiet",
    fwd_prefilt,
    rev_prefilt
  )
  
  if (!all(file.exists(c(cutadapt_data_barcode$trimmed_path))) && !barcode_params$already_trimmed) {
    cutadapt_output <- furrr::future_map(command_args, ~ system2(cutadapt, args = .x))
  }
  
  else if (!all(file.exists(c(cutadapt_data_barcode$trimmed_path))) && barcode_params$already_trimmed) {
    cutadapt_output <- furrr::future_map(command_args, ~ system2(cutadapt, args = .x))
    
    for (i in seq_along(fwd_untrim)) {
      fwd_untrim_reads <- ShortRead::readFastq(fwd_untrim[i])
      rev_untrim_reads <- ShortRead::readFastq(rev_untrim[i])
      ShortRead::writeFastq(fwd_untrim_reads, fwd_trim[i], mode = 'a')
      cat("Already trimmed forward reads were appended to trimmed read directory, and they are located here:", fwd_trim[i], "\n")
      ShortRead::writeFastq(rev_untrim_reads, rev_trim[i], mode = 'a')
      cat("Already trimmed reverse reads were appended to trimmed read directory, and they are located here:", rev_trim[i], "\n")
    }
  }
}

#' Wrapper function for plotQualityProfile function
#'
#' @inheritParams dada2::plotQualityProfile
#'
#' @param cutadapt_data FASTQ read files trimmed of primers
#' @param output_directory_path The path to the directory where resulting files
#'   are output
#' @param seed Set seed for reproducibility
#'
#' @return Dada2 wrapper function for making quality profiles for each sample
#'
#' @keywords internal
plot_qc <- function(cutadapt_data, output_directory_path, n = 500000, barcode_params) {
  #just retrieve all plots for first sample
  for (i in unique(cutadapt_data$sample_name)) {
    name1 = paste0('readqual_pretrim_plot_', i, '.pdf')
    pdf_path = file.path(output_directory_path, name1)
    set.seed = barcode_params$seed
    if (!file.exists(pdf_path)) {
      sample_info = cutadapt_data$trimmed_path[cutadapt_data$sample_name == i]
      quality_plots <- dada2::plotQualityProfile(sample_info, n)
      ggplot2::ggsave(
        quality_plots,
        filename = name1,
        path = output_directory_path,
        width = 8,
        height = 8
      )
    }
  }
}

#' Wrapper function for filterAndTrim function from DADA2, to be used after
#' primer trimming
#'
#' @inheritParams dada2::filterAndTrim
#' @param output_directory_path The path to the directory where resulting files
#'   are output
#' @param cutadapt_data_barcode Metabarcode-specific FASTQ read files trimmed of primers
#'
#' @return Filtered and trimmed reads
#'
#' @keywords internal
filter_and_trim <- function(output_directory_path,
                            temp_directory_path,
                            cutadapt_data_barcode,
                            barcode_params,
                            barcode) {
  filtered_read_dir <- file.path(temp_directory_path, "filtered_sequences")
  cutadapt_data_barcode$filtered_path <-
    file.path(
      filtered_read_dir,
      paste0(
        cutadapt_data_barcode$file_id,
        "_",
        cutadapt_data_barcode$primer_name,
        ".fastq.gz"
      )
    )
  if (!all(file.exists(cutadapt_data_barcode$filtered_path))) {
    seed = barcode_params$seed
    filter_results <-
      dada2::filterAndTrim(
        fwd = cutadapt_data_barcode$trimmed_path[cutadapt_data_barcode$direction == "Forward"],
        filt = cutadapt_data_barcode$filtered_path[cutadapt_data_barcode$direction == "Forward"],
        rev = cutadapt_data_barcode$trimmed_path[cutadapt_data_barcode$direction == "Reverse"],
        filt.rev = cutadapt_data_barcode$filtered_path[cutadapt_data_barcode$direction == "Reverse"],
        maxN = barcode_params$maxN,
        maxEE = c(barcode_params$maxEE_forward, barcode_params$maxEE_reverse),
        truncLen = c(barcode_params$truncLen_forward, barcode_params$truncLen_reverse),
        truncQ = barcode_params$truncQ,
        minLen = barcode_params$minLen,
        maxLen = barcode_params$maxLen,
        minQ = barcode_params$minQ,
        trimLeft = barcode_params$trimLeft,
        trimRight = barcode_params$trimRight,
        rm.phix = TRUE,
        compress = TRUE,
        matchIDs = TRUE,
        multithread = barcode_params$multithread,
        verbose = barcode_params$verbose,
        OMP = TRUE,
        n = 1e+05,
        id.sep = "\\s",
        rm.lowcomplex = barcode_params$rm.lowcomplex,
        orient.fwd = NULL,
        qualityType = "Auto",
        id.field = NULL
      )
    filter_results <- dplyr::as_tibble(filter_results)
    filter_results_path <-
      file.path(temp_directory_path, paste0("Filter_results_", barcode, ".RData"))
    save(filter_results, file = filter_results_path)
  }
}

#' Get primer counts for reach sample after primer removal and trimming steps
#'
#' @param primer_data The primer data frame created in orient_primers function
#' @param cutadapt_data FASTQ read files trimmed of primers
#' @param output_directory_path The path to the directory where resulting files
#'   are output
#'
#' @return Table of read counts across each sample
#'
#' @keywords internal
get_post_trim_hits <- function(primer_data, cutadapt_data, output_directory_path) {
  post_trim_hit_data <- tidyr::gather(
    primer_data,
    key = "orientation",
    value = "sequence",
    forward,
    f_compt,
    f_rev,
    f_rc,
    reverse,
    r_compt,
    r_rev,
    r_rc
  )
  
  # Function to count the number of reads in which the primer is found
  post_trim_hits <- function(primer, path) {
    nhits <- Biostrings::vcountPattern(primer, ShortRead::sread(ShortRead::readFastq(path), pattern = function(x) x), fixed = FALSE)
    return(sum(nhits > 0))
  }
  
  post_primer_hit_data_csv_path <- file.path(output_directory_path, "primer_hit_data_posttrim.csv")
  post_primer_hit_counts <- furrr::future_map(cutadapt_data$filtered_path,
                                              function(a_path)
                                                purrr::map_dbl(post_trim_hit_data$sequence, post_trim_hits, path = a_path))
  
  names(post_primer_hit_counts) <- paste0(cutadapt_data$file_id, "_", cutadapt_data$primer_name)
  post_trim_hit_data <- dplyr::bind_cols(post_trim_hit_data, dplyr::as_tibble(post_primer_hit_counts))
  readr::write_csv(post_trim_hit_data, post_primer_hit_data_csv_path)
  
  
  make_posttrim_primer_plot <- function(post_trim_hits,
                                        output_directory_path) {
    post_trim_hits <- post_trim_hits[-(3)]
    post_trim_hits$primer_type <- paste(post_trim_hits$primer_name, post_trim_hits$orientation)
    new_primer_hits <- post_trim_hits[-(1:2)]
    new_primer_hits <- new_primer_hits[, c("primer_type", setdiff(colnames(new_primer_hits), "primer_type"))]
    only_counts <- new_primer_hits[-(1)]
    colnames(only_counts) <- paste("Col", colnames(only_counts), sep = "-")
    total_nums <- apply(only_counts, 1, function(row) sum(as.numeric(row), na.rm = TRUE))
    new_primer_hits$Total <- paste(total_nums)
    total_primers <- new_primer_hits[, c("primer_type", "Total")]
    total_primers$Total <- as.numeric(total_primers$Total)
    plot <- ggplot2::ggplot(data = total_primers, ggplot2::aes(x = primer_type, y = Total)) +
      ggplot2::geom_bar(stat = "identity", width = 0.8, fill = "seagreen3") +
      ggplot2::geom_text(ggplot2::aes(label = Total), vjust = -0.5, color = "black", size = 3) +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "Number of primers found by barcode and orientation", x = "Primer Type", y = "Total")+
      ggplot2::theme_minimal()
    
    print(plot)
    ggplot2::ggsave(plot, filename = file.path(output_directory_path, "posttrim_primer_plot.pdf"), width = 8, height = 8)
    return(invisible(plot))
  }
  
  make_posttrim_primer_plot(post_trim_hit_data, output_directory_path)
}

#' Wrapper script for plotQualityProfile after trim steps and primer removal.
#'
#' @inheritParams dada2::plotQualityProfile
#' @param cutadapt_data FASTQ read files trimmed of primers
#' @param output_directory_path The path to the directory where resulting files
#'   are output
#'
#' @return Quality profiles of reads after primer trimming
#'
#' @keywords internal
plot_post_trim_qc <- function(cutadapt_data, output_directory_path, n = 500000, barcode_params) {
  for (i in unique(cutadapt_data$sample_name)) {
    name2 = paste0('readqual_posttrim_plot_', i, '.pdf')
    pdf_path = file.path(output_directory_path, name2)
    if (!file.exists(pdf_path)) {
      sample_info2 = cutadapt_data$filtered_path[cutadapt_data$sample_name == i]
      set.seed = barcode_params$seed
      quality_plots2 <- dada2::plotQualityProfile(sample_info2, n)
      ggplot2::ggsave(
        quality_plots2,
        filename = name2,
        path = output_directory_path,
        width = 8,
        height = 8
      )
    }
  }
}

#' Main command to trim primers using Cutadapt and core DADA2 functions   
#' 
#' @importFrom utils modifyList read.table stack
#' 
#' @param analysis_setup An object containing directory paths and data tables,
#'   produced by the `prepare_reads` function
#' @param cutadapt_path Path to the Cutadapt program.
#' @param overwrite_existing Logical, indicating whether to remove or overwrite
#'   existing files and directories from previous runs. Default is `FALSE`.
#'   
#' @return Trimmed reads, primer counts, quality plots, and ASV matrix.
#' 
#' @details If samples are comprised of two different barcodes (like ITS1 and rps10), reads will also be demultiplexed prior to DADA2 trimming steps. 
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Remove remaining primers from raw reads, demultiplex pooled barcoded samples, 
#' # and then trim reads based on specific DADA2 parameters
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
#' }
cut_trim <- function(analysis_setup,
                     cutadapt_path,
                     overwrite_existing = FALSE) {
  
  data_tables <- analysis_setup$data_tables
  output_directory_path <- analysis_setup$directory_paths$output_directory
  temp_directory_path <- analysis_setup$directory_paths$temp_directory
  
  patterns_to_check <- c(
    "primer_hit_data_posttrim.csv", 
    "posttrim_primer_plot.pdf",
    "readqual*"
  )
  
  files_exist <- sapply(patterns_to_check, function(pattern) {
    full_pattern <- file.path(output_directory_path, pattern)
    any(file.exists(list.files(path = output_directory_path, pattern = pattern, full.names = TRUE, recursive = TRUE)))
  })
  
  if (any(files_exist) && !overwrite_existing) {
    message("Existing data detected: Primer counts and N's may have been removed from previous runs. Loading existing output. To perform a new analysis, specify overwrite_existing = TRUE.")
    return(invisible())
  } else if (overwrite_existing) {
    warning("Existing analysis files found. Overwriting existing files.")
    
    # Remove existing files and directories
    patterns_to_remove <- c(
      "primer_hit_data_posttrim.csv", 
      "posttrim_primer_plot.pdf",
      "readqual*"
    )
    
    for (pattern in patterns_to_remove) {
      full_pattern <- file.path(output_directory_path, pattern)
      files_to_remove <- list.files(path = output_directory_path, pattern = pattern, full.names = TRUE)
      
      if (length(files_to_remove) > 0) {
        file.remove(files_to_remove)
      }
    }
    
    patterns_to_remove_temp <- c(
      "Filter_results_*"
    )
    
    for (pattern in patterns_to_remove_temp) {
      full_pattern <- file.path(temp_directory_path, pattern)
      files_to_remove <- list.files(path = temp_directory_path, pattern = pattern, full.names = TRUE)
      
      if (length(files_to_remove) > 0) {
        file.remove(files_to_remove)
      }
    }
    
    temp_untrimmed <- file.path(temp_directory_path, "untrimmed_sequenced")
    temp_trimmed <- file.path(temp_directory_path, "trimmed_sequences")
    temp_filtered <- file.path(temp_directory_path, "filtered_sequences")
    
    if (dir.exists(temp_untrimmed)) {
      unlink(list.files(temp_untrimmed, full.names = TRUE), recursive = TRUE)
    }
    if (dir.exists(temp_trimmed)) {
      unlink(list.files(temp_trimmed, full.names = TRUE), recursive = TRUE)
    }
    if (dir.exists(temp_filtered)) {
      unlink(list.files(temp_filtered, full.names = TRUE), recursive = TRUE)
    }
    
    subdirectory_names <- c("filtered_sequences", "trimmed_sequences", "untrimmed_sequences")
    
    for (seqdir_name in subdirectory_names) {
      seqdir_path <- file.path(temp_directory_path, seqdir_name)
      
      if (dir.exists(seqdir_path)) {
        unlink(list.files(seqdir_path, full.names = TRUE), recursive = TRUE)
      }
    }
  }
  
  default_params <- list(
    minCutadaptlength=0,
    maxEE_forward = Inf,
    maxEE_reverse = Inf,
    truncQ = 2,
    minLen = 20,
    maxLen = Inf,
    truncLen_forward = 0,
    truncLen_reverse = 0,
    maxN = 0,
    minQ = 0,
    trimLeft = 0,
    trimRight = 0,
    rm.phix = TRUE,
    multithread = FALSE,
    verbose = FALSE,
    qualityType = "Auto",
    OMP = TRUE,
    n = 1e+05,
    id.sep = "\\s",
    rm.lowcomplex = 0,
    orient.fwd = NULL,
    id.field = NULL,
    seed = NULL,
    overwrite_existing = FALSE
  )
  
  unique_primers <- unique(data_tables$primer_data$primer_name)
  
  for (barcode in unique_primers) {
    barcode_params <- dplyr::filter(data_tables$parameters, primer_name == barcode)
    
    if (nrow(barcode_params) > 0) {
      barcode_params <- as.list(barcode_params)
      barcode_params <- utils::modifyList(default_params, barcode_params)
      
      cutadapt_data_barcode <- subset(data_tables$cutadapt_data, primer_name == barcode)
      if (nrow(cutadapt_data_barcode) > 0 && !all(file.exists(c(cutadapt_data_barcode$trimmed_path)))) {
        run_cutadapt(
          cutadapt_path,
          cutadapt_data_barcode,
          barcode_params,
          minCutadaptlength = barcode_params$minCutadaptlength
        )
        
        if (length(barcode_params) > 0) {
          barcode_params <- as.list(barcode_params)
          filter_and_trim(
            output_directory_path,
            temp_directory_path,
            cutadapt_data_barcode,
            barcode_params, 
            barcode
          )
        }
        
        # Call plot_post_trim_qc within the loop
        quality_plots <- plot_qc(cutadapt_data_barcode, output_directory_path, barcode_params = barcode_params)
        post_primer_hit_data <- get_post_trim_hits(data_tables$primer_data, cutadapt_data_barcode, output_directory_path)
        plot_post_trim_qc(cutadapt_data_barcode, output_directory_path, barcode_params = barcode_params)
      }
    }
  }
}