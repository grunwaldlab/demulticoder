#' Prepare reads for primer trimming using Cutadapt
#' @param data_directory User-specified directory path where the user has placed
#'   raw FASTQ (forward and reverse reads), metadata.csv, and
#'   primerinfo_params.csv files. Default is "data".
#' @param output_directory User-specified directory for outputs. Default is
#'   "output".
#' @param tempdir_path Path to a temporary directory. If `NULL`, a temporary
#'   directory path will be identified using the `tempdir()` command.
#' @param tempdir_id ID for temporary directories. Default is
#'   "demulticoder_run". The user can provide any helpful ID, whether it be a
#'   date or specific name for the run.
#' @param multithread Logical, indicating whether to use multithreading for
#'   certain operations. Default is `FALSE`.
#' @param overwrite_existing Logical, indicating whether to remove or overwrite
#'   existing files and directories from previous runs.
#'
#' @return A list containing data tables, including metadata, primer sequences
#'   to search for based on orientation, paths for trimming reads, and
#'   user-defined parameters for all subsequent steps.
#'
#' @export
#'
#' @examples
#' # Pre-filter raw reads and parse metadata and primer_information to prepare 
#' # for primer trimming and filter
#' prepare_reads(
#'   data_directory = system.file("extdata", package = "demulticoder"), 
#'   output_directory = tempdir(),
#'   tempdir_path = tempdir(),
#'   tempdir_id = "demulticoder_run_temp", 
#'   overwrite_existing = FALSE
#' )
prepare_reads <- function(data_directory = "data",
                          output_directory = "output",
                          tempdir_path = NULL,
                          tempdir_id = "demulticoder_run",
                          multithread = FALSE,
                          overwrite_existing = FALSE) {
  directory_paths <-
    setup_directories(data_directory, output_directory, tempdir_path, tempdir_id)
  output_directory_path <- directory_paths$output_directory
  data_directory_path <- directory_paths$data_directory
  temp_directory_path <- directory_paths$temp_directory
  primers_params_path <- directory_paths$primers_params_path
  metadata_path <- directory_paths$metadata_path
  
  existing_files <- list.files(output_directory_path)
  
  if (!overwrite_existing && length(existing_files) > 0) {
    existing_analysis_table <-
      file.path(temp_directory_path, "analysis_setup_obj.RData")
    if (file.exists(existing_analysis_table)) {
      load(existing_analysis_table)
      message(
        "Existing data detected: Primer counts and N's may have been removed from previous runs. Loading existing output. To perform a new analysis, specify overwrite_existing = TRUE."
      )
    } else {
      warning(
        "Existing analysis setup tables are not found. The 'prepare_reads' function was rerun"
      )
    }
  } else {
    if (dir.exists(output_directory_path)) {
      unlink(output_directory_path, recursive = TRUE)
    }
    if (dir.exists(temp_directory_path)) {
      unlink(temp_directory_path, recursive = TRUE)
    }
  }
  
  if (!dir.exists(output_directory_path)) {
    dir.create(output_directory_path, recursive = TRUE)
  }
  
  if (!dir.exists(temp_directory_path)) {
    dir.create(temp_directory_path, recursive = TRUE)
  }
  
  primer_data <- orient_primers(primers_params_path)
  parameters <- read_parameters(primers_params_path)
  metadata <- prepare_metadata_table(metadata_path, primer_data)
  fastq_data <-
    read_prefilt_fastq(data_directory_path, multithread, temp_directory_path)
  pre_primer_hit_data <-
    get_pre_primer_hits(primer_data, fastq_data, output_directory_path)
  cutadapt_data <-
    make_cutadapt_tibble(fastq_data, metadata, temp_directory_path)
  data_tables <- list(
    cutadapt_data = cutadapt_data,
    primer_data = primer_data,
    fastq_data = fastq_data,
    parameters = parameters,
    metadata = metadata
  )
  
  analysis_setup <-
    list(data_tables = data_tables, directory_paths = directory_paths)
  analysis_setup_path <-
    file.path(temp_directory_path,
              paste0("analysis_setup_obj", ".RData"))
  save(analysis_setup, file = analysis_setup_path)
  return(analysis_setup)
}
#' Set up directory paths for subsequent analyses
#'
#' This function sets up the directory paths for subsequent analyses. It checks
#' whether the specified output directories exist or creates them if they don't.
#' The function also provides paths to primer and metadata files within the data
#' directory.
#'
#' @param data_directory User-specified directory path where the user has placed
#'   raw FASTQ (forward and reverse reads), metadata.csv, and
#'   primerinfo_params.csv files. Default is "data".
#' @param output_directory User-specified directory for outputs. Default is
#'   "output".
#' @param tempdir_path Path to a temporary directory. If `NULL`, a temporary
#'   directory path will be identified using the `tempdir()` command.
#' @param tempdir_id ID for temporary directories. Default is
#'   "demulticoder_run". The user can provide any helpful ID, whether it be a
#'   date or specific name for the run.
#'
#' @return A list with paths for data, output, temporary directories, primer,
#'   and metadata files.
#'
#' @keywords internal
setup_directories <- function(data_directory = "data",
                              output_directory = "output",
                              tempdir_path = NULL,
                              tempdir_id = "demulticoder_run") {
  primers_params_path <-
    file.path(data_directory, "primerinfo_params.csv")
  metadata_file_path <- file.path(data_directory, "metadata.csv")
  
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }
  
  output_path <- file.path(output_directory)
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  if (is.null(tempdir_path)) {
    temp_path <- file.path(tempdir(), paste0(tempdir_id))
  } else {
    temp_path <- file.path(tempdir_path, paste0(tempdir_id))
  }
  
  if (!dir.exists(temp_path)) {
    dir.create(temp_path, recursive = TRUE)
  }
  
  outputs <- (
    list(
      data_directory = data_directory,
      output_directory = output_directory,
      temp_directory = temp_path,
      primers_params_path = primers_params_path,
      metadata_path = metadata_file_path
    )
  )
}

#' Read metadata file from user and combine and reformat it, given primer data.
#' Included in a larger function prepare_reads.
#'
#' @param primer_data A data frame of oriented primer information returned from
#'   the orient_primers function.
#' @param metadata_path The path to the metadata file.
#'
#' @return A dataframe containing the merged metadata and primer data.
#'
#' @keywords internal
prepare_metadata_table <-
  function(metadata_file_path, primer_data) {
    metadata <- readr::read_csv(metadata_file_path)
    metadata <- merge(metadata, primer_data, by = "primer_name")
    metadata <- metadata[order(metadata$sample_name), ]
    if ("primer_name" %in% colnames(metadata)) {
      new_metadata_cols <-
        c("primer_name", colnames(metadata)[colnames(metadata) != "primer_name"])
      metadata <- metadata[new_metadata_cols]
      metadata$samplename_barcode <-
        paste0(metadata$sample_name, "_", metadata$primer_name)
      barcode_col <- match("samplename_barcode", colnames(metadata))
      metadata <-
        metadata[, c(barcode_col, seq_along(metadata)[-barcode_col])]
    } else {
      stop("Please make sure that there is a 'primer_name' column in your metadata table.")
    }
    return(metadata)
  }

#' Take in user's forward and reverse sequences and creates the complement,
#' reverse, reverse complement of primers in one data frame
#'
#' @param primers_params_path A path to the CSV file that holds the primer
#'   information.
#'
#' @return A data frame with oriented primer information.
#'
#' @keywords internal
orient_primers <- function(primers_params_path) {
  primer_data_path <- file.path(primers_params_path)
  primer_data <- readr::read_csv(primer_data_path)
  
  forward_primers <- primer_data[, c(1:2)]
  reverse_primers <- primer_data[, c(1, 3)]
  
  
  forward_primers$f_compt <-
    sapply(forward_primers$forward, function(x)
      toString(Biostrings::complement(Biostrings::DNAString(x))))
  forward_primers$f_rev <-
    sapply(forward_primers$forward, function(x)
      toString(Biostrings::reverse(Biostrings::DNAString(x))))
  forward_primers$f_rc <-
    sapply(forward_primers$forward, function(x)
      toString(Biostrings::reverseComplement(Biostrings::DNAString(x))))
  
  reverse_primers$r_compt <-
    sapply(reverse_primers$reverse, function(x)
      toString(Biostrings::complement(Biostrings::DNAString(x))))
  reverse_primers$r_rev <-
    sapply(reverse_primers$reverse, function(x)
      toString(Biostrings::reverse(Biostrings::DNAString(x))))
  reverse_primers$r_rc <-
    sapply(reverse_primers$reverse, function(x)
      toString(Biostrings::reverseComplement(Biostrings::DNAString(x))))
  
  #add back together
  primer_data <-
    merge(forward_primers, reverse_primers, by = "primer_name")
  return(primer_data)
}

#' Take in user's DADA2 parameters and make a dataframe for downstream steps
#'
#' @param primers_params_path A path to the CSV file that holds the primer
#'   information.
#'
#' @return A data frame with information on the DADA2 parameters.
#'
#' @keywords internal
read_parameters <- function(primers_params_path) {
  params_data_path <- file.path(primers_params_path)
  param_data <- readr::read_csv(params_data_path)
  parameters <- param_data[, c(1, 4:ncol(param_data))]
  return(parameters)
}

#' Takes in the FASTQ files from the user and creates a data frame with the
#' paths to files that will be created and used in the future. Included in a
#' larger 'read_prefilt_fastq' function.
#'
#' @param data_directory_path The path to the directory containing the FASTQ,
#'   metadata, and primer_info files
#' @param temp_directory_path User-defined temporary directory to place reads
#'   throughout the workflow.
#'
#' @return A data frame with the FASTQ file paths, primer orientations and
#'   sequences, and parsed sample names
#'
#' @keywords internal
read_fastq <- function(data_directory_path,
                       temp_directory_path) {
  fastq_paths <-
    list.files(data_directory_path, pattern = "\\.fastq(|\\.gz)$")
  
  file_id <- sub("\\.fastq(|\\.gz)$", "", fastq_paths)
  sample_name = gsub(fastq_paths, pattern = "_R1|_R2\\.fastq|.fastq\\.gz|.gz$", replacement = "")
  
  direction <-
    ifelse(grepl("_R1", fastq_paths), "Forward", "Reverse")
  
  fastq_data_directory_paths <-
    file.path(data_directory_path, fastq_paths)
  temp_data_path <- file.path(temp_directory_path, fastq_paths)
  
  fastq_data <-
    data.frame(file_id,
               sample_name,
               direction,
               fastq_data_directory_paths,
               temp_data_path)
  
  return(fastq_data)
}

#' Matching Order Primer Check
#'
#' @param fastq_data A data frame with the FASTQ file paths, the direction of
#'   the sequences, and names of sequences
#'
#' @return None
#'
#' @keywords internal
primer_check <- function(fastq_data) {
  paired_file_paths <-
    fastq_data[fastq_data$sample_name == fastq_data$sample_name[1], "fastq_data_directory_paths"]
  
  get_read_names <- function(path) {
    seqs <- read.table(path, sep = "\n", stringsAsFactors = FALSE)
    sub(".+$", "", seqs$V1)
  }
  
  stopifnot(all(
    get_read_names(paired_file_paths[1]) == get_read_names(paired_file_paths[2])
    
  ))
}


#' A function for calling read_fastq, primer_check, and remove_ns functions.
#' This will process and edit the FASTQ and make them ready for the trimming of
#' primers with Cutadapt. Part of a larger 'prepare_reads' function.
#'
#' @inheritParams dada2::filterAndTrim
#'
#' @param data_directory_path The path to the directory containing the FASTQ,
#'   metadata.csv, and primerinfo_params.csv files
#' @param temp_directory_path User-defined temporary directory to output
#'   unfiltered, trimmed, and filtered read directories throughout the workflow
#'
#' @return Returns filtered reads that have no Ns
#'
#' @keywords internal
read_prefilt_fastq <-
  function(data_directory_path = data_directory_path,
           multithread = FALSE,
           temp_directory_path) {
    fastq_data <- read_fastq(data_directory_path, temp_directory_path)
    primer_check(fastq_data)
    fastq_data <-
      remove_ns(fastq_data, multithread = multithread, temp_directory_path)
    return(fastq_data)
  }

#' Wrapper function for core DADA2 filter and trim function for first filtering
#' step
#'
#' @inheritParams dada2::filterAndTrim
#'
#' @param fastq_data A data frame with the fastq file paths, the direction of
#'   the sequences, and names of sequences metadata, and primer_info files
#'
#' @param metadata A metadata containing the concatenated metadata and primer
#'   data
#' @param temp_directory_path User-defined temporary directory to output
#'   unfiltered, trimmed, and filtered read directories throughout the workflow
#' @return Return prefiltered reads with no Ns
#'
#' @keywords internal
remove_ns <-
  function(fastq_data,
           multithread = TRUE,
           temp_directory_path) {
    prefiltered_read_dir <-
      file.path(temp_directory_path, "prefiltered_sequences")
    fastq_data$prefiltered_path <-
      file.path(prefiltered_read_dir,
                base::basename(fastq_data$temp_data_path))
    
    if (!all(file.exists(fastq_data$prefiltered_path))) {
      dada2::filterAndTrim(
        fwd = fastq_data[fastq_data$direction == "Forward",][["fastq_data_directory_paths"]],
        filt = fastq_data[fastq_data$direction == "Forward",][["prefiltered_path"]],
        rev = fastq_data[fastq_data$direction == "Reverse",][["fastq_data_directory_paths"]],
        filt.rev = fastq_data[fastq_data$direction == "Reverse",][["prefiltered_path"]],
        multithread = multithread,
        maxN = 0
      )
    }
    return(fastq_data)
  }

#' Get primer counts for reach sample before primer removal and trimming steps
#'
#' @param output_directory_path The path to the directory where resulting files
#'   are output
#' @param primer_data The primer data data frame created in orient_primers
#'   function
#' @param fastq_data A data frame with the FASTQ file paths, the direction of
#'   the sequences, and names of sequences
#'
#' @return The number of reads in which the primer is found
#'
#' @keywords internal
get_pre_primer_hits <-
  function(primer_data,
           fastq_data,
           output_directory_path) {
    primer_hit_data <-
      tidyr::gather(
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
    
    primer_hits <- function(primer, path) {
      fastq_data <- ShortRead::readFastq(path)
      nhits <- Biostrings::vcountPattern(primer, ShortRead::sread(fastq_data), fixed = FALSE)
      return(sum(nhits > 0))
    }
    
    primer_hit_data_csv_path <-
      file.path(output_directory_path, "primer_hit_data_pretrim.csv")
    
    if (file.exists(primer_hit_data_csv_path)) {
      primer_hit_data <- readr::read_csv(primer_hit_data_csv_path)
      
    } else {
      primer_hit_counts <- furrr::future_map(fastq_data$prefiltered_path,
                                      function (a_path)
                                        purrr::map_dbl(primer_hit_data$sequence, primer_hits, path = a_path))
      
      names(primer_hit_counts) <- paste0(fastq_data$file_id)
      
      primer_hit_data <-
        dplyr::bind_cols(primer_hit_data, dplyr::as_tibble(primer_hit_counts))
      readr::write_csv(primer_hit_data, primer_hit_data_csv_path)
    }
    
    make_primer_hit_plot <- function(primer_hits,
                                     output_directory_path) {
      primer_hits <- primer_hits[-(3)]
      primer_hits$primer_type <-
        paste(primer_hits$primer_name, primer_hits$orientation)
      new_primer_hits <- primer_hits[-(1:2)]
      new_primer_hits <-
        new_primer_hits[, c("primer_type", setdiff(colnames(new_primer_hits), "primer_type"))]
      only_counts <- new_primer_hits[-(1)]
      colnames(only_counts) <-
        paste("Col", colnames(only_counts), sep = "-")
      total_nums <-
        apply(only_counts, 1, function(row)
          sum(as.numeric(row), na.rm = TRUE))
      new_primer_hits$Total <- paste(total_nums)
      total_primers <- new_primer_hits[, c("primer_type", "Total")]
      total_primers$Total <- as.numeric(total_primers$Total)
      plot <-
        ggplot2::ggplot(data = total_primers, ggplot2::aes(x = primer_type, y = Total)) +
        ggplot2::geom_bar(stat = "identity",
                 width = 0.8,
                 fill = "seagreen3") +
        ggplot2::geom_text(
          ggplot2::aes(label = Total),
          vjust = -0.5,
          color = "black",
          size = 3
        ) +
        ggplot2::coord_flip() +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Number of primers found by barcode and orientation", x = "Primer Type", y = "Total")
      
      print(plot)
      ggplot2::ggsave(
        plot,
        filename = file.path(output_directory_path, "pretrim_primer_plot.pdf"),
        width = 8,
        height = 8
      )
      return(invisible(plot))
    }
    
    make_primer_hit_plot(primer_hit_data, output_directory_path)
  }

#' Prepare for primmer trimming with Cutaapt. Make new sub-directories and
#' specify paths for the trimmed and untrimmed reads
#'
#' @param fastq_data A path to FASTQ files for analysis metadata, and
#'   primer_info files
#' @param metadata Loaded metadata pairing the user's metadata file with the
#'   primer data
#' @param temp_directory_path User-defined temporary directory to output
#'   unfiltered, trimmed, and filtered read directories throughout the workflow
#'
#' @return Returns a larger data frame containing paths to temporary read
#'   directories, which is used as input when running Cutadapt
#'
#' @keywords internal
make_cutadapt_tibble <-
  function(fastq_data,
           metadata,
           temp_directory_path) {
    cutadapt_data <- merge(metadata, fastq_data, by = "sample_name")
    
    trimmed_read_dir <-
      file.path(temp_directory_path, "trimmed_sequences")
    if (!dir.exists(trimmed_read_dir)) {
      dir.create(trimmed_read_dir)
    }
    
    cutadapt_data$trimmed_path <-
      file.path(
        trimmed_read_dir,
        paste0(
          cutadapt_data$file_id,
          "_",
          cutadapt_data$primer_name,
          ".fastq.gz"
        )
      )
    
    untrimmed_read_dir <-
      file.path(temp_directory_path, "untrimmed_sequences")
    if (!dir.exists(untrimmed_read_dir)) {
      dir.create(untrimmed_read_dir)
    }
    cutadapt_data$untrimmed_path <-
      file.path(
        untrimmed_read_dir,
        paste0(
          cutadapt_data$file_id,
          "_",
          cutadapt_data$primer_name,
          ".fastq.gz"
        )
      )
    
    filtered_read_dir <-
      file.path(temp_directory_path, "filtered_sequences")
    if (!dir.exists(filtered_read_dir)) {
      dir.create(filtered_read_dir)
    }
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
    return(cutadapt_data)
  }
