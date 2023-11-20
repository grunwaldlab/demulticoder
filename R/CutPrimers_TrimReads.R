#' Main command to trim primers based on Cutadapt and DADA2 functions. If samples contain pooled barcodes, reads will also be demultiplexed. 
#'
#' @param analysis_setup A list containing directory paths and data tables, produced by the `prepare_reads` function.
#' @param cutadapt_path A path to the Cutadapt program.
#' @param rawSeqTab_fileName A filename as which the raw sequence table will be saved.
#' @param abundMatrix_fileName A filename as which the abundance matrix will be saved.
#' @param overwrite_existing Logical, indicating whether to remove or overwrite existing files and directories from previous runs. If set to TRUE, specific output files 
#' @return Reads trimmed of primers and filtered, primer counts after running Cutadapt, quality plots after poor quality reads are trimmed or removed, and the ASV matrix.
#' @export
#' @inheritParams plot_qc
#' @inheritParams filter_and_trim
#' @inheritParams run_cutadapt
#' 
#' @examples
#' analysis_setup <- prepare_reads(
#'   data_directory = system.file("extdata", package = "your_package_name"),
#'   output_directory = tempdir(),
#'   tempdir_id = "run1",
#'   overwrite_existing = FALSE)
#'
#' # Main function to trim primers based on Cutadapt and DADA2 functions
#' cut_trim(
#'   analysis_setup,
#'   cutadapt_path = "/opt/homebrew/bin/cutadapt",
#'   overwrite_existing = FALSE)


cut_trim <- function(analysis_setup,
                     cutadapt_path,
                     maxEE = Inf,
                     truncQ = 2,
                     minLen = 20,
                     maxLen = Inf,
                     truncLen = 0,
                     maxN = 0,
                     minQ = 0,
                     trimLeft = 0,
                     trimRight = 0,
                     rm.phix = TRUE,
                     multithread = FALSE,
                     matchIDs = FALSE,
                     verbose = FALSE,
                     qualityType = "Auto",
                     OMP = TRUE,
                     n = 1e+05,
                     id.sep = "\\s",
                     rm.lowcomplex = 0,
                     orient.fwd = NULL,
                     id.field = NULL,
                     minCutadaptlength = 50,
                     overwrite_existing = FALSE) {
  
  dir_paths <- analysis_setup$dir_paths
  data_tables <- analysis_setup$data_tables
  directory_path <- dir_paths$output_directory
  data_path <- dir_paths$data_directory
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
    
    patterns_to_remove_temp <- "filter_results.RData"
    
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
    
  run_cutadapt(cutadapt_path,
               data_tables$cutadapt_data,
               minCutadaptlength = minCutadaptlength)
  quality_plots <- plot_qc(data_tables$cutadapt_data, directory_path)
  filter_results <-
    filter_and_trim(
      directory_path,
      directory_path_temp,
      data_tables$cutadapt_data,
      maxEE = maxEE,
      truncQ = truncQ,
      minLen = minLen,
      maxLen = maxLen,
      multithread = multithread,
      matchIDs = matchIDs,
      verbose = verbose,
      qualityType = qualityType,
      OMP = OMP,
      n = n,
      id.sep = id.sep,
      trimLeft = trimLeft,
      trimRight = trimRight,
      rm.lowcomplex = rm.lowcomplex,
      orient.fwd = orient.fwd,
      id.field = id.field
    )
  post_primer_hit_data <-
    get_post_primer_hits(data_tables$primer_data,
                         data_tables$cutadapt_data,
                         directory_path)
  post_primer_plot <-
    make_primer_hit_plot(
      post_primer_hit_data,
      data_tables$fastq_data,
      directory_path,
      "posttrim_primer_plot.pdf"
    )
  quality_plots2 <-
    plot_post_trim_qc(data_tables$cutadapt_data, directory_path)
}


#' Core function for running cutadapt
#' Core function for running cutadapt
#'
#' @param cutadapt_path A path to the cutadapt program.
#' @param cutadapt_data Directory_data folder with trimmed and filtered reads for each sample.
#' @param minCutadaptlength Read lengths that are lower than this threshold will be discarded. Default is 50.
#' @return Trimmed read.
#' @keywords internal
run_cutadapt <- function(cutadapt_path,
                         cutadapt_data,
                         minCutadaptlength = 20) {
    cutadapt <- cutadapt_path
    tryCatch(
      system2(cutadapt, args = "--version"),
      warning = function(w) {
        stop("cutadapt cannot be found on PATH. Check if it is installed?")
      }
    )
    #simplify
    #may need to demultiplex first-meed to assess this
    cutadapt <- path.expand(cutadapt_path)
    minCutadaptlength = minCutadaptlength
    R1_flags = unique(paste("-g", cutadapt_data$forward, "-a", cutadapt_data$r_rc))
    R2_flags = unique(paste("-G", cutadapt_data$reverse, "-A", cutadapt_data$f_rc))
    fwd_trim = cutadapt_data[cutadapt_data$direction == "Forward",][["trimmed_path"]]
    rev_trim = cutadapt_data[cutadapt_data$direction == "Reverse",][["trimmed_path"]]
    fwd_untrim = cutadapt_data[cutadapt_data$direction == "Forward",][["untrimmed_path"]]
    rev_untrim = cutadapt_data[cutadapt_data$direction == "Reverse",][["untrimmed_path"]]
    fwd_prefilt = cutadapt_data[cutadapt_data$direction == "Forward",][["prefiltered_path"]]
    rev_prefilt = cutadapt_data[cutadapt_data$direction == "Reverse",][["prefiltered_path"]]
    #consider parameters
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
    if (!all(file.exists(c(cutadapt_data$trimmed_path)))) {
      cutadapt_output <-
        furrr::future_map(command_args, ~ system2(cutadapt, args = .x))
    }
  }

#' Wrapper function for plotQualityProfile function
#'
#' @inheritParams dada2::plotQualityProfile
#' @param cutadapt_data directory_data folder with trimmed and filtered reads for each sample
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primer_info files
#' @inheritParams plotQualityProfile
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
#' metadata, and primer_info files
#' @param cutadapt_data directory_data folder with trimmed and filtered reads for each sample
#' @return Filtered and trimmed reads
#' @keywords internal
#'
filter_and_trim <-
  function(directory_path,
           directory_path_temp,
           cutadapt_data,
           maxEE = Inf,
           truncQ = 2,
           minLen = 20,
           maxLen = Inf,
           truncLen = 0,
           maxN = 0,
           minQ = 0,
           trimLeft = 0,
           trimRight = 0,
           rm.phix = TRUE,
           multithread = FALSE,
           matchIDs = FALSE,
           verbose = FALSE,
           qualityType = "Auto",
           OMP = TRUE,
           n = 1e+05,
           id.sep = "\\s",
           rm.lowcomplex = 0,
           orient.fwd = NULL,
           id.field = NULL) {
    filtered_read_dir <- file.path(directory_path_temp, "filtered_sequences")
    cutadapt_data$filtered_path <-
      file.path(
        filtered_read_dir,
        paste0(
          cutadapt_data$file_id,
          "_",
          cutadapt_data$primer_name,
          ".fastq.gz"
        )
      )
    if (!all(file.exists(cutadapt_data$filtered_path))) {
      filter_results <-
        dada2::filterAndTrim(
          fwd = cutadapt_data$trimmed_path[cutadapt_data$direction == "Forward"],
          filt = cutadapt_data$filtered_path[cutadapt_data$direction == "Forward"],
          rev =  cutadapt_data$trimmed_path[cutadapt_data$direction == "Reverse"],
          filt.rev = cutadapt_data$filtered_path[cutadapt_data$direction == "Reverse"],
          maxN = maxN,
          maxEE = c(maxEE, maxEE),
          truncLen = c(truncLen, truncLen),
          truncQ = truncQ,
          minLen = minLen,
          maxLen = maxLen,
          minQ = minQ,
          trimLeft = 0,
          trimRight = 0,
          rm.phix = rm.phix,
          compress = TRUE,
          matchIDs = FALSE,
          multithread = multithread,
          verbose = verbose,
          OMP = TRUE,
          n = 1e+05,
          id.sep = "\\s",
          rm.lowcomplex = 0,
          orient.fwd = NULL,
          qualityType = "Auto",
          id.field = NULL
        )
      filter_results <- as_tibble(filter_results)
      filter_results_path <-
        file.path(directory_path_temp, "filter_results.RData")
      save(filter_results, file = filter_results_path)
      print(colMeans(filter_results))
    }
  }

#' Get primer counts fo reach sample after primer removal and trimming steps
#'
#' @param primer_data The primer data tibble created in orient_primers function
#' @param cutadapt_data directory_data folder with trimmed and filtered reads for each sample
#' @return Table of read counts across each sample
#' @keywords internal
#'

get_post_primer_hits <-
  function(primer_data,
           cutadapt_data,
           directory_path) {
    post_primer_hit_data <-
      gather(
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
    
    #from DADA2
    post_primer_hits <- function(primer, path) {
      # Counts number of reads in which the primer is found
      nhits <-
        vcountPattern(primer, sread(readFastq(path)), fixed = FALSE)
      return(sum(nhits > 0))
    }
    
    post_primer_hit_data_csv_path <-
      file.path(directory_path, "primer_hit_data_posttrim.csv")
    #if a file exists in there, then write the path to it
    if (file.exists(post_primer_hit_data_csv_path)) {
      post_primer_hit_data <- read_csv(post_primer_hit_data_csv_path)
      #if the file doesn't exist in primer hit data
      #map applied a function to each part of that specific vector
    } else {
      post_primer_hit_counts <- future_map(cutadapt_data$filtered_path,
                                           function (a_path)
                                             map_dbl(post_primer_hit_data$sequence, post_primer_hits, path = a_path))
      #gets or sets the name of an object
      names(post_primer_hit_counts) <-
        paste0(cutadapt_data$file_id, "_", cutadapt_data$primer_name)
      #primer hit data will be a tibble with the columns of primer hit data plus the primer hit counts
      #with names pulled from the fastq_data file
      post_primer_hit_data <-
        bind_cols(post_primer_hit_data, as_tibble(post_primer_hit_counts))
      print(post_primer_hit_data)
      write_csv(post_primer_hit_data, post_primer_hit_data_csv_path)
    }
    return(post_primer_hit_data)
  }

#' Wrapper script for plotQualityProfile after trim steps and primer removal.
#'
#' @inheritParams dada2::plotQualityProfile
#' @param cutadapt_data directory_data folder with trimmed and filtered reads for each sample
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primer_info files
#' @return Quality profiles of reads after primer trimming
#' @keywords internal
plot_post_trim_qc <-
  lot_post_trim_qc <-
  function(cutadapt_data, directory_path, n = 500000) {
    #just retrieve all plots for first sample
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