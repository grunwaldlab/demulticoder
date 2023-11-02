#' Set up directory paths for subsequent analyses
#'
#' This function sets up the paths for the analysis.
#' It checks whether the specified output directories exist or creates them if they don't. 
#' The function also provides paths to primer and metadata files within the data directory.
#'
#' @param data_directory Directory for data files. Default is "data".
#' @param output_directory Directory for analysis outputs. Default is "outputs".
#' @param tempdir_id ID for temporary directories, combined with the current date. Default is "run1".
#' @return A list with paths for data, output, temporary directories, primer, and metadata files.
#' @keywords internal
#' 
setup_directories <- function(data_directory = "data", 
                              output_directory = "outputs", 
                              tempdir_id = "run1") {
  
  data_primer_path <- file.path(data_directory, "primer_info.csv")
  data_metadata_path <- file.path(data_directory, "metadata.csv")
  
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }
  
  output_path <- file.path(output_directory)
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  temp_path <- file.path(tempdir(), paste0(tempdir_id, "_", Sys.Date()))
  if (!dir.exists(temp_path)) {
    dir.create(temp_path, recursive = TRUE)
  }
  
  outputs <- (list(data_directory = data_directory, 
              output_directory = output_directory, 
              temp_directory = temp_path, 
              primer_path = data_primer_path, 
              metadata_path = data_metadata_path))
}

#' Prepare reads for primer trimming using Cutadapt
#'
#' This function prepares sequencing reads for primer trimming using Cutadapt by setting up necessary directories, reading in metadata and primer data, and performing pre-trimming quality control steps.
#'
#' @inheritParams dada2::filterAndTrim
#' @param overwrite_existing Logical, indicating whether to remove or overwrite existing files and directories from previous runs. If set to TRUE, specific output files 
#' @param data_directory Directory for data files. Default is "data".
#' @param output_directory Directory for outputs. Default is "output".
#' @param tempdir_id ID for temporary directories, combined with the current date. Default is "run1".
#' @return A list with data tables for trimming and directory paths.
#' @export 
#' @examples
#' outputs <- prepare_reads(maxN = 0, data_directory = "data", output_directory = "output")

prepare_reads <-
  function(data_directory = "data", 
           output_directory = "output", 
           tempdir_id = "run1",
           maxN = 0, 
           multithread = FALSE, 
           overwrite_existing = FALSE) {
    
    dir_paths <- setup_directories(data_directory, output_directory, tempdir_id)
    directory_path <- dir_paths$output_directory
    data_path <- dir_paths$data_directory
    directory_path_temp <- dir_paths$temp_directory
    primer_path <- dir_paths$primer_path
    metadata_path <- dir_paths$metadata_path
    
    if (overwrite_existing) {
      if (dir.exists(directory_path)) {
        unlink(directory_path, recursive = TRUE)
      }
      if (dir.exists(directory_path_temp)) {
        unlink(directory_path_temp, recursive = TRUE)
      }
      
    }
    
    if (!dir.exists(directory_path)) {
      dir.create(directory_path, recursive = TRUE)
    }
    
    if (!dir.exists(directory_path_temp)) {
      dir.create(directory_path_temp, recursive = TRUE)
    }
    
    primer_data <- orient_primers(primer_path)
    metadata <- prepare_metadata_table(metadata_path, primer_data)
    fastq_data <- read_prefilt_fastq(data_path, maxN, multithread, directory_path_temp)
    pre_primer_hit_data <-
      get_pre_primer_hits(primer_data, fastq_data, directory_path)
    pre_primer_plot <-
      make_primer_hit_plot(pre_primer_hit_data,
                           fastq_data,
                           directory_path,
                           "pretrim_primer_plot.pdf")
    cutadapt_data <-
      make_cutadapt_tibble(fastq_data, metadata, directory_path_temp)
    data_tables <-
      list(
        cutadapt_data = cutadapt_data,
        primer_data = primer_data,
        fastq_data = fastq_data,
        metadata = metadata
      )
    analysis_setup <-list(data_tables = data_tables, dir_paths = dir_paths)
    assign("analysis_setup", analysis_setup, envir = .GlobalEnv)
    return(analysis_setup)
  }

#' Read in the metadata from user and combine it with the primer data.
#' Included in a larger function prepare_reads.
#' @param primer_path a path to the csv file that holds the primer
#' information
#' @param metadata_path The path to the metadata file
#' @return A dataframe containing the merged metadata and primer
#' data
#' @keywords internal
#'
prepare_metadata_table <- function(metadata_path, primer_data) {
  metadata <- read_csv(metadata_path)
  metadata <- merge(metadata, primer_data, by = "primer_name")
  metadata <- metadata[order(metadata$sample_name), ]
  if ("primer_name" %in% colnames(metadata)) {
    new_metadata_cols <- c("primer_name", colnames(metadata)[colnames(metadata) != "primer_name"])
    metadata <- metadata[new_metadata_cols]
    metadata$samplename_barcode <- paste0(metadata$sample_name, "_", metadata$primer_name)
    barcode_col <- match("samplename_barcode", colnames(metadata))
    metadata <- metadata[, c(barcode_col, seq_along(metadata)[-barcode_col])]
  } else {
    stop("Please make sure that there is a 'primer_name' column in your metadata table.")
  }
  return(metadata)
}

#' Take in user's primers and creates the complement, reverse,
#' reverse complement of primers in one tibble
#' @param primer_path a path to the csv file that holds the primer
#' information
#' @return A data frame with oriented primer information
#' @keywords internal
orient_primers <- function(primer_path) {
  primer_data_path <- file.path(primer_path)
  primer_data <- read_csv(primer_data_path)
  
  forward_primers <- primer_data[, c(1:2)]
  reverse_primers <- primer_data[, c(1, 3)]
  
  forward_primers$f_compt <- sapply(forward_primers$forward, function(x)
    toString(Biostrings::complement(DNAString(x))))
  forward_primers$f_rev <- sapply(forward_primers$forward, function(x)
    toString(Biostrings::reverse(DNAString(x))))
  forward_primers$f_rc <- sapply(forward_primers$forward, function(x)
    toString(Biostrings::reverseComplement(DNAString(x))))
  
  reverse_primers$r_compt <- sapply(reverse_primers$reverse, function(x)
    toString(Biostrings::complement(DNAString(x))))
  reverse_primers$r_rev <- sapply(reverse_primers$reverse, function(x)
    toString(Biostrings::reverse(DNAString(x))))
  reverse_primers$r_rc <- sapply(reverse_primers$reverse, function(x)
    toString(Biostrings::reverseComplement(DNAString(x))))
  
  #add back together
  primer_data <- merge(forward_primers, reverse_primers, by = "primer_name")
  
  return(primer_data)
}

#' Takes in the fastq files from the user and creates a tibble with
#' the paths to files that will be created and used in the future.
#' Included in a larger 'read_prefilt_fastq' function
#'
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primer_info files
#' @param directory_path_temp User-defined temporary directory to place reads throughout the workflow. 
#' @return A tibble with the fastq file paths, the direction of the
#' @keywords internal
#' sequences, and names of sequences

read_fastq <- function(data_path, directory_path_temp) {
  
  fastq_paths <- list.files(data_path, pattern = "\\.fastq(|\\.gz)$")

  file_id <- sub("\\.fastq(|\\.gz)$", "", fastq_paths)
  sample_name = gsub(fastq_paths, pattern = "_R1|_R2\\.fastq|.fastq\\.gz|.gz$", replacement = "")
  
  direction <- ifelse(grepl("_R1", fastq_paths), "Forward", "Reverse")

  directory_data_path <- file.path(data_path, fastq_paths)
  temp_data_path <- file.path(directory_path_temp, fastq_paths)
  
  fastq_data <- data.frame(file_id, sample_name, direction, directory_data_path, temp_data_path)
  
  return(fastq_data)
}

#' Matching Order Primer Check
#' @param fastq_data A tibble with the fastq file paths, the direction of
#' the sequences, and names of sequences
#' @return None
#' @keywords internal

primer_check <- function(fastq_data) {
  paired_file_paths <- fastq_data[fastq_data$sample_name == fastq_data$sample_name[1], "directory_data_path"]
  
  get_read_names <- function(path) {
    seqs <- read.table(path, sep = "\n", stringsAsFactors = FALSE)
    sub(".+$", "", seqs$V1)
  }
  
  stopifnot(all(
    get_read_names(paired_file_paths[1]) == get_read_names(paired_file_paths[2])
    
  ))
}


#' A function for calling read_fastq, primer_check, and remove_ns functions. This will process and edit the FASTQ and make them ready for the trimming of primers with Cutadapt. Part of a larger 'prepare_reads' function.
#'
#' @inheritParams dada2::filterAndTrim
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primer_info files
#' @param directory_path_temp User-defined temporary directory to place reads throughout the workflow
#' @return Returns filtered reads that have no Ns
#' @keywords internal

read_prefilt_fastq <- function(data_path, maxN = 0, multithread = FALSE, directory_path_temp) {
  fastq_data <- read_fastq(data_path, directory_path_temp)
  primer_check(fastq_data)
  fastq_data <- remove_ns(fastq_data, maxN, multithread = multithread, directory_path_temp)
  return(fastq_data)
}

#' Wrapper function for core DADA2 filter and trim function for first filtering step
#'
#' @inheritParams dada2::filterAndTrim
#' @param fastq_data A tibble with the fastq file paths, the direction of the sequences, and names of sequences
#' metadata, and primer_info files
#' @param metadata A metadata containing the concatenated metadata and primer data
#' @param directory_path_temp User-defined temporary directory to place reads throughout the workflow
#' @return Return prefiltered reads with no Ns
#' @keywords internal

remove_ns <-
  function(fastq_data,
           maxN = 0,
           multithread = TRUE, 
           directory_path_temp) {
    prefiltered_read_dir <-
      file.path(directory_path_temp, "prefiltered_sequences")
    fastq_data$prefiltered_path <-
      file.path(prefiltered_read_dir,
                base::basename(fastq_data$temp_data_path))
    
    if (!all(file.exists(fastq_data$prefiltered_path))) {
      dada2::filterAndTrim(
        fwd = fastq_data[fastq_data$direction == "Forward",][["directory_data_path"]],
        filt = fastq_data[fastq_data$direction == "Forward",][["prefiltered_path"]],
        rev = fastq_data[fastq_data$direction == "Reverse",][["directory_data_path"]],
        filt.rev = fastq_data[fastq_data$direction == "Reverse",][["prefiltered_path"]],
        maxN = maxN,
        multithread = multithread
      )
    }
    return(fastq_data)
  }

#' Get primer counts for reach sample before primer removal and trimming steps
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primer_info files
#' @param primer_data The primer data tibble created in orient_primers function
#' @param fastq_data A tibble with the fastq file paths, the direction of the sequences, and names of sequences
#' @return The number of reads in which the primer is found
#' @keywords internal

get_pre_primer_hits <-
  function(primer_data, fastq_data, directory_path) {
    primer_hit_data <-
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
    
    primer_hits <- function(primer, path) {
      nhits <-
        vcountPattern(primer, sread(readFastq(path)), fixed = FALSE)
      return(sum(nhits > 0))
    }
    
    primer_hit_data_csv_path <-
      file.path(directory_path, "primer_hit_data_pretrim.csv")

    if (file.exists(primer_hit_data_csv_path)) {
      primer_hit_data <- read_csv(primer_hit_data_csv_path)
      print(primer_hit_data)

    } else {
      primer_hit_counts <- future_map(fastq_data$prefiltered_path,
                                      function (a_path)
                                        map_dbl(primer_hit_data$sequence, primer_hits, path = a_path))
     
      names(primer_hit_counts) <- paste0(fastq_data$file_id)

      primer_hit_data <-
        bind_cols(primer_hit_data, as_tibble(primer_hit_counts))
      print(primer_hit_data)
      write_csv(primer_hit_data, primer_hit_data_csv_path)
    }
    return(primer_hit_data)
    
  }

#' Make a barplot of primers identified on reads
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primer_info files
#' @param primer_hits A number of reads in which the primer is found
#' @param fastq_data A tibble with the fastq file paths, the direction of
#' the sequences, and names of sequences
#' @param plot_name A filename under which a PDF file of the plot will be saved as
#' @return Returns a barplot with read counts
#' @keywords internal

make_primer_hit_plot <- function(primer_hits,
                                 fastq_data,
                                 directory_path,
                                 plot_name) {
  
  primer_hits <- primer_hits[-(3)]
  primer_hits$primer_type <- paste(primer_hits$primer_name, primer_hits$orientation)
  new_primer_hits <- primer_hits[-(1:2)]
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
    labs(title = "Number of primers found by barcode and orientation", x = "Primer Type", y = "Total")
  
  print(plot)
  
  ggsave(plot, filename = file.path(directory_path, paste0(plot_name)), width = 8, height = 8)
  
  return(invisible(plot))
}

#' Prepare for primmer trimming with Cutaapt. Make new sub-directories
#' and specify paths for the trimmed and untrimmed reads
#' @param fastq_data A path to FASTQ files for analysis
#' metadata, and primer_info files
#' @param metadata A metadata containing the concatenated metadata and primer data
#' @param directory_path_temp User-defined temporary directory to place reads throughout the workflow.
#' @return Returns a tibble that is used as input when running Cutadapt
#' @keywords internal
make_cutadapt_tibble <-
  function(fastq_data, metadata, directory_path_temp) {
    cutadapt_data <- merge(metadata, fastq_data, by = "sample_name")
    
    trimmed_read_dir <- file.path(directory_path_temp, "trimmed_sequences")
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
      file.path(directory_path_temp, "untrimmed_sequences")
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
      file.path(directory_path_temp, "filtered_sequences")
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
