#Prepare reads. A wrapper function to prepare reads for trimming using Cutadapt. Counts of primers on reads will be output.
#' Main command prepare reads for primer trimming
#'
#' @param directory_path The path to the directory containing the fastq,
#' @param directory_path_temp User-defined temporary directory to place reads throughout the workflow
#' metadata, and primer_info files
#' @param primer_path a path to the csv file that holds the primer
#' information
#' @param metadata_path The path to the metadata file
#' @inheritParams read_prefilt_fastq
#' @return A list containinpre_primer_hit_data the following data tables:
#' \itemize{
#'   \item \code{cutadapt_data}: The cutadapt data table.
#'   \item \code{primer_data}: The primer data table.
#'   \item \code{fastq_data}: The fastq data table.
#'   \item \code{metadata}: The metadata table.
#' }
#' @export prepare_reads
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
#' maxN = 0,
#' )
#'

prepare_reads <-
  function(directory_path,
           directory_path_temp,
           primer_path,
           metadata_path,
           maxN = 0,
           multithread = FALSE) {
    primer_data <- orient_primers(primer_path)
    metadata <- prepare_metadata_table(metadata_path, primer_data)
    fastq_data <- read_prefilt_fastq(directory_path, maxN, multithread, directory_path_temp)
    pre_primer_hit_data <-
      get_pre_primer_hits(primer_data, fastq_data, directory_path)
    pre_primer_plot <-
      make_primer_hit_plot(pre_primer_hit_data,
                           fastq_data,
                           directory_path,
                           "pre_primer_plot.pdf")
    cutadapt_data <-
      make_cutadapt_tibble(fastq_data, metadata, directory_path_temp)
    data_tables <-
      list(
        cutadapt_data = cutadapt_data,
        primer_data = primer_data,
        fastq_data = fastq_data,
        metadata = metadata
      )
    return(data_tables)
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
  
  #separate forward and reverse to make various primers
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

read_fastq <- function(directory_path, directory_path_temp) {
  
  fastq_paths <- list.files(directory_path, pattern = "\\.fastq(|\\.gz)$")
  
  # Extract file_id and sample_name using regex
  file_id <- sub("\\.fastq(|\\.gz)$", "", fastq_paths)
  sample_name = gsub(fastq_paths, pattern = "_R1|_R2\\.fastq|.fastq\\.gz|.gz$", replacement = "")
  
  # Determine direction based on file name
  direction <- ifelse(grepl("_R1", fastq_paths), "Forward", "Reverse")

  directory_data_path <- file.path(directory_path, fastq_paths)
  temp_data_path <- file.path(directory_path_temp, fastq_paths)
  
  # Create the data frame
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

read_prefilt_fastq <- function(directory_path, maxN = 0, multithread = FALSE, directory_path_temp) {
  fastq_data <- read_fastq(directory_path, directory_path_temp)
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
    
    #from DADA2
    primer_hits <- function(primer, path) {
      # Counts number of reads in which the primer is found
      nhits <-
        vcountPattern(primer, sread(readFastq(path)), fixed = FALSE)
      return(sum(nhits > 0))
    }
    
    primer_hit_data_csv_path <-
      file.path(directory_path, "primer_hit_data_pre_trim.csv")
    #if a file exists in there, then write the path to it
    if (file.exists(primer_hit_data_csv_path)) {
      primer_hit_data <- read_csv(primer_hit_data_csv_path)
      print(primer_hit_data)
      #if the file doesn't exist in primer hit data
      #map applied a function to each part of that specific vector
    } else {
      primer_hit_counts <- future_map(fastq_data$prefiltered_path,
                                      function (a_path)
                                        map_dbl(primer_hit_data$sequence, primer_hits, path = a_path))
      #gets or sets the name of an object
      names(primer_hit_counts) <- paste0(fastq_data$file_id)
      #primer hit data will be a tibble with the columns of primer hit data plus the primer hit counts
      #with names pulled from the fastq_data file
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
  
  # Removing the sequence in the primer_hits tibble
  primer_hits <- primer_hits[-(3)]
  
  # Concatenating the two beginning columns
  primer_hits$primer_type <- paste(primer_hits$primer_name, primer_hits$orientation)
  
  # Subsetting just the concatenated column
  new_primer_hits <- primer_hits[-(1:2)]
  
  # Move name to the first column
  new_primer_hits <- new_primer_hits[, c("primer_type", setdiff(colnames(new_primer_hits), "primer_type"))]
  
  # Easier to work without the character names
  only_counts <- new_primer_hits[-(1)]
  
  # Add a prefix to the columns so they are easier to mutate
  colnames(only_counts) <- paste("Col", colnames(only_counts), sep = "-")
  
  # Add them all up
  total_nums <- apply(only_counts, 1, function(row) sum(as.numeric(row), na.rm = TRUE))
  
  # Add back to new_primer_hits
  new_primer_hits$Total <- paste(total_nums)
  
  # Subset just the names and the totals
  total_primers <- new_primer_hits[, c("primer_type", "Total")]
  
  # Convert the total column to numeric
  total_primers$Total <- as.numeric(total_primers$Total)
  
  # Create a bar chart using ggplot2
  plot <- ggplot(data = total_primers, aes(x = primer_type, y = Total)) +
    geom_bar(stat = "identity", width = 0.8, fill = "seagreen3") +
    geom_text(aes(label = Total), vjust = -0.5, color = "black", size = 3) +
    coord_flip() +
    labs(title = "Number of primers found by barcode and orientation", x = "Primer Type", y = "Total")
  
  # Print the plot
  print(plot)
  
  # Save the plot to a file
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
    # merge metadata and fastq_data by 'sample_name'
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
