#' Main command to trim primers based on Cutadapt and DADA2 functions. If samples contain pooled barcodes, reads will also be demultiplexed. 
#'
#' @param analysis_setup A list containing directory paths and data tables, produced by the `prepare_reads` function.
#' @param cutadapt_path A path to the Cutadapt program.
#' @param overwrite_existing Logical, indicating whether to remove or overwrite existing files and directories from previous runs.
#' @return Reads trimmed of primers and filtered, primer counts after running Cutadapt, quality plots after poor quality reads are trimmed or removed, and the ASV matrix.
#' @export
#' @examples
#' Remove remaining primers from raw reads, demultiplex pooled barcoded samples, and then trim reads based on specific DADA2 parameters
#' prepare_reads(maxN = 0, data_directory = "~/demulticoder/inst/extdata", output_directory = "~/testing_package", tempdir_id = "run1", overwrite_existing = TRUE)
#' cut_trim(analysis_setup,cutadapt_path="/opt/homebrew/bin/cutadapt", overwrite_existing = TRUE)
cut_trim <- function(analysis_setup,
                     cutadapt_path,
                     overwrite_existing = FALSE) {
  
  dir_paths <- analysis_setup$dir_paths
  data_tables <- analysis_setup$data_tables
  directory_path <- dir_paths$output_directory
  directory_path_temp <- dir_paths$temp_directory
  
  if (overwrite_existing) {
    patterns_to_remove <- c(
      "primer_hit_data_posttrim.csv", 
      "posttrim_primer_plot.pdf",
      "readqual*"
    )
    
    for (pattern in patterns_to_remove) {
      full_pattern <- file.path(directory_path, pattern)
      files_to_remove <- list.files(path = directory_path, pattern = pattern, full.names = TRUE)
      
      if (length(files_to_remove) > 0) {
        file.remove(files_to_remove)
      }
    }
    
    patterns_to_remove_temp <- "Filter_results.RData"
    
    for (pattern in patterns_to_remove_temp) {
      full_pattern <- file.path(directory_path_temp, pattern)
      files_to_remove <- list.files(path = directory_path_temp, pattern = pattern, full.names = TRUE)
      
      if (length(files_to_remove) > 0) {
        file.remove(files_to_remove)
      }
    }
    
    subdirectory_names <- c("filtered_sequences", "trimmed_sequences", "untrimmed_sequences")
    
    for (seqdir_name in subdirectory_names) {
      seqdir_path <- file.path(directory_path_temp, seqdir_name)
      
      if (dir.exists(seqdir_path)) {
        files_to_remove <- list.files(seqdir_path, full.names = TRUE)
        file.remove(files_to_remove)
      }
    }
  }
  
  default_params <- list(
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
    overwrite_existing = FALSE
  )
  
  unique_primers <- unique(data_tables$primer_data$primer_name)
  
  for (barcode in unique_primers) {
    barcode_params <- filter(data_tables$parameters, primer_name == barcode)
    
    if (nrow(barcode_params) > 0) {
      barcode_params <- as.list(barcode_params)
      
      # Merge default_params with barcode_params, favoring barcode-specific ones
      barcode_params <- modifyList(default_params, barcode_params)
      
      cutadapt_data_barcode <- subset(data_tables$cutadapt_data, primer_name == barcode)
      if (nrow(cutadapt_data_barcode) > 0 && !all(file.exists(c(cutadapt_data_barcode$trimmed_path)))) {
        run_cutadapt(
          cutadapt_path,
          cutadapt_data_barcode,
          barcode_params,
          minCutadaptlength = barcode_params$minCutadaptlength)
        
        quality_plots <- plot_qc(cutadapt_data_barcode, directory_path)
        
        if (length(barcode_params) > 0) {
          barcode_params <- as.list(barcode_params)
          filter_and_trim(
            directory_path,
            directory_path_temp,
            cutadapt_data_barcode,
            barcode_params, 
            barcode)
        }
      }
    }
  }
  post_primer_hit_data <- get_post_trim_hits(data_tables$primer_data, data_tables$cutadapt_data, directory_path)
  quality_plots2 <- plot_post_trim_qc(data_tables$cutadapt_data, directory_path)
}

#' Core function for running cutadapt

#' @param cutadapt_path A path to the cutadapt program.
#' @param cutadapt_data Directory_data folder with trimmed and filtered reads for each sample.
#' @param minCutadaptlength Read lengths that are lower than this threshold will be discarded. Default is 50.
#' @return Trimmed read.
#' @keywords internal
#' Core function for running cutadapt
#' Core function for running cutadapt
#'
run_cutadapt <- function(cutadapt_path,
                         cutadapt_data_barcode,
                         barcode_params,
                         minCutadaptlength) {
  
  cutadapt <- cutadapt_path
  tryCatch(
    system2(cutadapt, args = "--version"),
    warning = function(w) {
      stop("cutadapt cannot be found on PATH. Check if it is installed?")
    }
  )
  
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
  
  # Print command arguments
  print("Command Arguments:")
  print(command_args)
  
  # Run cutadapt if trimmed files don't exist
  if (!all(file.exists(c(cutadapt_data_barcode$trimmed_path))) && !barcode_params$already_trimmed) {
    cutadapt_output <- furrr::future_map(command_args, ~ system2(cutadapt, args = .x))
    # Print Cutadapt output
    print("Cutadapt Output:")
    print(cutadapt_output)
  }
  
  else if (!all(file.exists(c(cutadapt_data_barcode$trimmed_path))) && barcode_params$already_trimmed) {
    cutadapt_output <- furrr::future_map(command_args, ~ system2(cutadapt, args = .x))
    print("cutadapt_output")
    
    for (i in seq_along(fwd_untrim)) {
      fwd_untrim_reads <- readFastq(fwd_untrim[i])
      rev_untrim_reads <- readFastq(rev_untrim[i])
      writeFastq(fwd_untrim_reads, fwd_trim[i], mode = 'a')
      writeFastq(rev_untrim_reads, rev_trim[i], mode = 'a')
      cat("Contents appended to trimmed files.\n")
    }
  }
}
#' Wrapper function for plotQualityProfile function
#'
#' @inheritParams dada2::plotQualityProfile
#' @param cutadapt_data directory_data folder with trimmed and filtered reads for each sample
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primers and params files
#' @return Dada2 wrapper function for making quality profiles for each sample
#' @keywords internal

plot_qc <- function(cutadapt_data, directory_path, n = 500000) {
  #just retrieve all plots for first sample
  for (i in unique(cutadapt_data$sample_name)) {
    name1 = paste0('readqual_pretrim_plot_', i, '.pdf')
    pdf_path = file.path(directory_path, name1)
    if (!file.exists(pdf_path)) {
      sample_info = cutadapt_data$trimmed_path[cutadapt_data$sample_name == i]
      quality_plots <- dada2::plotQualityProfile(sample_info, n)
      ggplot2::ggsave(
        quality_plots,
        filename = name1,
        path = directory_path,
        width = 8,
        height = 8
      )
    }
  }
}

#' Wrapper function for filterAndTrim function from DADA2 after primer removal
#'
#' @inheritParams dada2::filterAndTrim
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primersinfo_params files
#' @param cutadapt_data_barcode directory_data folder with trimmed and filtered reads for each sample
#' @return Filtered and trimmed reads
#' @keywords internal
#'

filter_and_trim <- function(directory_path,
                            directory_path_temp,
                            cutadapt_data_barcode,
                            barcode_params,
                            barcode) {
  filtered_read_dir <- file.path(directory_path_temp, "filtered_sequences")
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
    filter_results <- as_tibble(filter_results)
    filter_results_path <-
      file.path(directory_path_temp, paste0("Filter_results_", barcode, ".RData"))
    save(filter_results, file = filter_results_path)
    print(colMeans(filter_results))
  }
}
#' Get primer counts for reach sample after primer removal and trimming steps
#'
#' @param primer_data The primer data tibble created in orient_primers function
#' @param cutadapt_data directory_data folder with trimmed and filtered reads for each sample
#' @param directory_path 
#' @return Table of read counts across each sample
#' @keywords internal
#'

get_post_trim_hits <- function(primer_data, cutadapt_data, directory_path) {
  post_trim_hit_data <- gather(
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
    nhits <- vcountPattern(primer, sread(readFastq(path)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  
  post_primer_hit_data_csv_path <- file.path(directory_path, "primer_hit_data_posttrim.csv")
  post_primer_hit_counts <- future_map(cutadapt_data$trimmed_path,
                                       function(a_path)
                                         map_dbl(post_trim_hit_data$sequence, post_trim_hits, path = a_path))
  
  names(post_primer_hit_counts) <- paste0(cutadapt_data$file_id, "_", cutadapt_data$primer_name)
  post_trim_hit_data <- bind_cols(post_trim_hit_data, as_tibble(post_primer_hit_counts))
  write_csv(post_trim_hit_data, post_primer_hit_data_csv_path)
  
  
  make_posttrim_primer_plot <- function(post_trim_hits,
                                        directory_path) {
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
    plot <- ggplot(data = total_primers, aes(x = primer_type, y = Total)) +
      geom_bar(stat = "identity", width = 0.8, fill = "seagreen3") +
      geom_text(aes(label = Total), vjust = -0.5, color = "black", size = 3) +
      coord_flip() +
      labs(title = "Number of primers found by barcode and orientation", x = "Primer Type", y = "Total")+
      theme_minimal()
    
    print(plot)
    ggsave(plot, filename = file.path(directory_path, "posttrim_primer_plot.pdf"), width = 8, height = 8)
    return(invisible(plot))
  }
  
  make_posttrim_primer_plot(post_trim_hit_data, directory_path)
}
# Assuming check_primer_counts is not needed anymore
# You can directly call the function with post_primer_hit_data


#' Get primer counts for reach sample after primer removal and trimming steps
#'
#' @param primer_data The primer data tibble created in orient_primers function
#' @param cutadapt_data directory_data folder with trimmed and filtered reads for each sample
#' @param directory_path 
#' @return Table of read counts across each sample
#' @keywords internal
#'

get_post_trim_hits <- function(primer_data, cutadapt_data, directory_path) {
  post_trim_hit_data <- gather(
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
    nhits <- vcountPattern(primer, sread(readFastq(path)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  
  post_primer_hit_data_csv_path <- file.path(directory_path, "primer_hit_data_posttrim.csv")
  post_primer_hit_counts <- future_map(cutadapt_data$trimmed_path,
                                       function(a_path)
                                         map_dbl(post_trim_hit_data$sequence, post_trim_hits, path = a_path))
  
  names(post_primer_hit_counts) <- paste0(cutadapt_data$file_id, "_", cutadapt_data$primer_name)
  post_trim_hit_data <- bind_cols(post_trim_hit_data, as_tibble(post_primer_hit_counts))
  write_csv(post_trim_hit_data, post_primer_hit_data_csv_path)
  
  
  make_posttrim_primer_plot <- function(post_trim_hits,
                                        directory_path) {
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
    plot <- ggplot(data = total_primers, aes(x = primer_type, y = Total)) +
      geom_bar(stat = "identity", width = 0.8, fill = "seagreen3") +
      geom_text(aes(label = Total), vjust = -0.5, color = "black", size = 3) +
      coord_flip() +
      labs(title = "Number of primers found by barcode and orientation", x = "Primer Type", y = "Total")+
      theme_minimal()
    
    print(plot)
    ggsave(plot, filename = file.path(directory_path, "posttrim_primer_plot.pdf"), width = 8, height = 8)
    return(invisible(plot))
  }
  
  make_posttrim_primer_plot(post_trim_hit_data, directory_path)
}
# Assuming check_primer_counts is not needed anymore
# You can directly call the function with post_primer_hit_data


#' Wrapper script for plotQualityProfile after trim steps and primer removal.
#'
#' @inheritParams dada2::plotQualityProfile
#' @param cutadapt_data directory_data folder with trimmed and filtered reads for each sample
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primersinfo_params files
#' @return Quality profiles of reads after primer trimming
#' @keywords internal
plot_post_trim_qc <- function(cutadapt_data, directory_path, n = 500000) {
  for (i in unique(cutadapt_data$sample_name)) {
    name2 = paste0('readqual_posttrim_plot_', i, '.pdf')
    pdf_path = file.path(directory_path, name2)
    if (!file.exists(pdf_path)) {
      sample_info2 = cutadapt_data$filtered_path[cutadapt_data$sample_name == i]
      quality_plots2 <- dada2::plotQualityProfile(sample_info2, n)
      ggplot2::ggsave(
        quality_plots2,
        filename = name2,
        path = directory_path,
        width = 8,
        height = 8
      )
    }
  }
}