#Prepare reads. A wrapper function to prepare reads for trimming using Cutadapt. Counts of primers on reads will be output.
#' Main command prepare reads for primer trimming
#'
#' @param directory_path A path to the intermediate folder and directory
#' @param primer_path The primer data tibble created
#' in orient_primers function
#' @param metadata_path A path to a metadata containing
#' the concatenated metadata and primer data
#' @inheritParams read_prefilt_fastq
#' @return A list containing data modified by cutadapt,
#' primer data, FASTQ data, and concatenated metadata and primer data
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
           primer_path,
           metadata_path,
           maxN = 0,
           multithread = FALSE,
           force = FALSE) {
    if(!force & file.exists("pre_primer_plot.pdf")) {
      unlink(c("pre_primer_plot.pdf", "primer_hit_data_pretrim.csv",
               "filtered_sequences", "prefiltered_sequences",
               "trimmed_sequences", "untrimmed_sequences"))
    }
    primer_data <- orient_primers(primer_path)
    metadata <- prepare_metadata_table(metadata_path, primer_data)
    fastq_data <-
      read_prefilt_fastq(directory_path, maxN = maxN, multithread = multithread)
    pre_primer_hit_data <-
      get_pre_primer_hits(primer_data, fastq_data, directory_path)
    pre_primer_plot <-
      make_primer_hit_plot(pre_primer_hit_data,
                           fastq_data,
                           directory_path,
                           "pre_primer_plot.pdf")
    cutadapt_data <-
      make_cutadapt_tibble(fastq_data, metadata, directory_path)
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
#' @param directory_path Folder path of read files and metadata file
#' @param metadata_path  A path to the metadata .csv file
#' @param primer_data The primer data tibble created in orient_primers
#' function
#' @return A metadata table containing the concatenated metadata and primer
#' data
#' @keywords internal
#'
prepare_metadata_table <- function(metadata_path, primer_data) {
  metadata <- read_csv(metadata_path) %>%
    left_join(primer_data, by = c("primer_name"))
  metadata <- metadata[order(metadata$sample_name), ]
  if ("primer_name" %in% colnames(metadata)) {
    # Check if primer_name is in metadata columns, if true, moves it to the first column
    new_metadata_cols <-
      c("primer_name", colnames(metadata)[colnames(metadata) != "primer_name"])
    metadata = metadata[new_metadata_cols]
    metadata$sample_nameBarcode <-
      paste0(metadata$sample_name, "_", metadata$primer_name)
    metadata <-
      metadata %>% relocate(sample_nameBarcode, .after = sample_name)
  }
  else {
    # Raises error if false
    stop("Please make sure that there is a 'primer_name' column in your metadata table.")
  }
  return(metadata)
}

#' Take in user's primers and creates the complement, reverse,
#' reverse complement of primers in one tibble
#' @param primer_path a path to the csv file that holds the primer
#' info for this project
#' @return A tibble that contains the forward and reverse primers
#' with all complements of primers
#' @keywords internal

orient_primers <- function(primer_path) {
  primer_data_path <- file.path(primer_path)
  primer_data <- read_csv(primer_data_path)

  #separate forward and reverse to make various primers
  forward_primers <- primer_data[, c(1:2)]
  toString(forward_primers)
  reverse_primers <- primer_data[, c(1, 3)]

  forward_primers$f_compt <-
    map_chr(forward_primers$forward, function(x)
      toString(Biostrings::complement(DNAString(x))))
  forward_primers$f_rev <-
    map_chr(forward_primers$forward, function(x)
      toString(Biostrings::reverse(DNAString(x))))
  forward_primers$f_rc <-
    map_chr(forward_primers$forward, function(x)
      toString(Biostrings::reverseComplement(DNAString(x))))

  reverse_primers$r_compt <-
    map_chr(reverse_primers$reverse, function(x)
      toString(Biostrings::complement(DNAString(x))))
  reverse_primers$r_rev <-
    map_chr(reverse_primers$reverse, function(x)
      toString(Biostrings::reverse(DNAString(x))))
  reverse_primers$r_rc <-
    map_chr(reverse_primers$reverse, function(x)
      toString(Biostrings::reverseComplement(DNAString(x))))

  #add back together
  primer_data <- forward_primers %>%
    left_join(reverse_primers, by = "primer_name")

  return(primer_data)
}

#' Takes in the fastq files from the user and creates a tibble with
#' the paths to files that will be created and used in the future.
#' Included in a larger 'read_prefilt_fastq' function
#'
#' @param directory_path Folder path of read files and metadata file
#' @return A tibble with the fastq file paths, the direction of the
#' @keywords internal
#' sequences, and names of sequences

read_fastq <- function(directory_path) {
  fastq_paths <- list.files(directory_path, pattern = "\\.fastq")
  #constructing a tibble - special data frame with improved behaviors.
  fastq_data <-
    tibble(
      file_id = sub(fastq_paths, pattern = "\\.fastq\\.gz$", replacement = ""),
      sample_name = gsub(fastq_paths, pattern = "_R1|_R2\\.fastq|.fastq\\.gz|.gz$", replacement = ""),
      #direction: the grep1 is looking to match the pattern, relates a TRUE or FALSE with it (aka 1 or 0) and adds 1 to find the
      #correct index to put Reverse or Forward
      direction = c("Reverse", "Forward")[BiocGenerics::grepl(fastq_paths, pattern = "_R1") + 1],
      directory_data_path = file.path(directory_path, fastq_paths)
    )

  return(fastq_data)
}

#' Matching Order Primer Check
#' @param fastq_data A tibble with the fastq file paths, the direction of
#' the sequences, and names of sequences
#' @return None
#' @keywords internal

primer_check <- function(fastq_data) {
  paired_file_paths <- fastq_data %>%
    filter(sample_name == gdata::first(sample_name)) %>%
    #pull is just like $ - you are accessing a variable inside the tibble
    pull(directory_data_path)

  #this is the creation of a function: needs a path passed to it
  #literally just checking if the files have the same number of ids/seq
  get_read_names <- function(path) {
    #from ShortRead: really convenient to just read in fasta
    seqs <- readFastq(path)
    #replaces all the seq ids as characters, then replaces everything in array with empty ""
    sub(as.character(seqs@id),
        pattern = ".+$",
        replacement = "")
  }
  #if the expression above is not true, it will return a error message
  #so if the forward and reverse reads are not in matching order
  #comparing the get read names at index 1 and 2 in the paired_file_paths
  #all is comparing the empty arrays and seeing if there is the same number in them
  stopifnot(all(
    get_read_names(paired_file_paths[1]) == get_read_names(paired_file_paths[2])
  ))
}

#' A function for calling read_fastq, primer_check, and remove_ns functions. This will process and edit the FASTQ and make them ready for the trimming of primers with Cutadapt. Part of a larger 'prepare_reads' function.
#'
#' @inheritParams dada2::filterAndTrim
#' @param directory_path Folder path of read files and metadata file
#' @param metadata Metadata containing the concatenated metadata and primer data
#' @return Returns filtered reads that have no Ns
#' @keywords internal

read_prefilt_fastq <-
  function(directory_path,
           maxN = 0,
           multithread = FALSE) {
    fastq_data <- read_fastq(directory_path)
    primer_check(fastq_data)
    fastq_data <-
      remove_ns(fastq_data, directory_path, maxN, multithread = multithread)
    return(fastq_data)
  }

#' Wrapper function for core DADA2 filter and trim function for first filtering step
#'
#' @inheritParams dada2::filterAndTrim
#' @param fastq_data A tibble with the fastq file paths, the direction of the sequences, and names of sequences
#' @param directory_path Folder path of read files and metadata file
#' @param metadata A metadata containing the concatenated metadata and primer data
#' @inheritParams filterAndTrim
#' @return Return prefiltered reads with no Ns
#' @keywords internal

remove_ns <-
  function(fastq_data,
           directory_path,
           maxN = 0,
           multithread = TRUE) {
    prefiltered_read_dir <-
      file.path(directory_path, "prefiltered_sequences")
    fastq_data$prefiltered_path <-
      file.path(prefiltered_read_dir,
                base::basename(fastq_data$directory_data_path))
    #if the files do not exist in the prefiltered path path (which clearly they don't)
    #raw path is the path to the actual fastq files in mock community
    #fwd takes from the raw path data
    #filt is the path to the output filtered files that we created
    #the sequences it chose to put into prefiltered path had to have no Ns in them
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

#' Get primer counts fo reach sample before primer removal and trimming steps
#' @param directory_path Folder path of read files and metadata file
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
      file.path(directory_path, "primer_hit_data_pretrim.csv")
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
#' @param directory_path Folder path of read files and metadata file
#' @param primer_hits A number of reads in which the primer is found
#' @param fastq_data A tibble with the fastq file paths, the direction of
#' the sequences, and names of sequences
#' @param plot_name A filename under which a PDF file of the plot will be saved as
#' @return Returns a barplot with read counts
#' @keywords internal

make_primer_hit_plot <-
  function(primer_hits,
           fastq_data,
           directory_path,
           plot_name) {
    #removing the sequence in the primer_hits tibble
    primer_hits <- primer_hits[-(3)]
    #concatonating the two beginning columns
    primer_hits$primer_type <-
      paste(primer_hits$primer_name, primer_hits$orientation)
    #subsetting just the concatenated column
    new_primer_hits <- primer_hits[-(1:2)]
    #move name to the first column
    new_primer_hits <- new_primer_hits %>%
      select(primer_type, everything())

    #easier to work without the character names
    only_counts <- new_primer_hits[-(1)]
    #add a prefix to the columns so they are easier to mutate
    colnames(only_counts) <-
      paste("Col", colnames(only_counts), sep = "-")

    #add them all up
    total_nums <- only_counts %>%
      rowwise() %>%
      mutate(Total = sum(across(starts_with("Col")), na.rm = TRUE))

    #add back to new_primer_hits
    new_primer_hits$Total <- paste(total_nums$Total)

    #subset just the names and the totals
    needed_cols <- c("primer_type", "Total")
    total_primers <- new_primer_hits[needed_cols]

    #convert the total column to numeric
    total_primers <-
      transform(total_primers, Total = as.numeric(Total))
    #a bar chart
    plot <- ggplot(data = total_primers, aes(x = primer_type, y = Total)) +
      geom_bar(stat = "identity",
               width = 1.0,
               fill = "seagreen3") +
      geom_text(aes(label = Total)) +
      coord_flip()
    print(plot)
    ggsave(
      plot,
      filename = plot_name,
      path = directory_path,
      width = 8,
      height = 8
    )
    return(plot)
  }

#' Prepare for primmer trimming with Cutaapt. Make new sub-directories
#' and specify paths for the trimmed and untrimmed reads
#'
#' @param fastq_data A path to FASTQ files for analysis
#' @param directory_path Folder path of read files and metadata file
#' @param metadata A metadata containing the concatenated metadata and primer data
#' @return Returns a tibble that is used as input when running Cutadapt
#' @keywords internal
make_cutadapt_tibble <-
  function(fastq_data, metadata, directory_path) {
    #new_fastq_data needed, why not just fastq_data
    cutadapt_data <- metadata %>%
      left_join(fastq_data, by = c("sample_name" = "sample_name"))
    trimmed_read_dir <- file.path(directory_path, "trimmed_sequences")
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
      file.path(directory_path, "untrimmed_sequences")
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
      file.path(directory_path, "filtered_sequences")
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
