#' Create intermediate path and directory
#'
#' @param working_dir_path A path to a directory containing reads and metadata/primer files
#'
#' @return A path to the intermediate folder and directory
#' @export
#'
#' @examples intermediate_path <- create_intermediate(directory_path)
create_intermediate <- function(working_dir_path){
  intermediate_read_dir <- file.path(working_dir_path,  "intermediate_data")
  if (! dir.exists(intermediate_read_dir)){
    dir.create(intermediate_read_dir)
  }
  intermediate_path <- file.path(intermediate_read_dir)
  return(intermediate_path)
}

#' Read in the metadata from user and combine it with the primer data. Included in a larger function prepare_reads. 
#'
#' @param metadata_path  A path to the metadata .csv file
#' @param primer_data The primer data tibble created in prepare_primers function
#'
#' @return A metadata containing the concatenated metadata and primer data
#' @export
#'
#' @examples metadata <- prepare_metadata(metadata_path, primer_data)
prepare_metadata <- function(metadata_path, primer_data){
  metadata <- read_csv(metadata_path) %>%
    left_join(primer_data, by = c("primer_name"))
  if ("primer_name" %in% colnames(metadata)) { # Check if primer_name is in metadata columns, if true, moves it to the first column
    new_metadata_cols <- c("primer_name", colnames(metadata)[colnames(metadata) != "primer_name"]) 
    metadata = metadata[new_metadata_cols]
  }
  else { # Raises error if false
    stop("Please make sure that there is a 'primer_name' column in your metadata table.")
  }
  #this could cause issues in other data sets, if the name of the sample is not the first column
  #names(metadata)[1] <- "sample_id"
  return(metadata)
}

#' Takes in the fastq files from the user and creates a tibble with the paths to files that will be created and used in the future. Included in a larger 'prepare_fastq' function
#'
#' @param raw_path The path to the fastq files from user
#'
#' @return A tibble with the fastq file paths, the direction of the sequences, and names of sequences
#' @export
#'
#' @examples fastq_data <- read_fastq(raw_path)
read_fastq <- function(raw_path){
  fastq_paths <- list.files(raw_path, pattern = "\\.fastq")
  #constructing a tibble - special data frame with improved behaviors.
  fastq_data <- tibble(file_id = sub(fastq_paths, pattern = "\\.fastq\\.gz$", replacement = ""),
                       sample_id = sub(fastq_paths, pattern = "_.+$", replacement = ""),
                       #direction: the grep1 is looking to match the pattern, relates a TRUE or FALSE with it (aka 1 or 0) and adds 1 to find the
                       #correct index to put Reverse or Forward
                       direction = c( "Reverse", "Forward")[BiocGenerics::grepl(fastq_paths, pattern = "_R1") + 1],
                       raw_data_path = file.path(raw_path, fastq_paths))

  return(fastq_data)
}


#' Matching Order Primer Check
#'
#' @param fastq_data A tibble with the fastq file paths, the direction of the sequences, and names of sequences
#'
#' @return None
#' @export
#'
#' @examples Part of a larger function 'prepare_fastq' function. 
matching_order_primer_check <- function(fastq_data){
  paired_file_paths <- fastq_data %>%
    filter(sample_id == gdata::first(sample_id)) %>%
    #pull is just like $ - you are accessing a variable inside the tibble
    pull(raw_data_path)

  #this is the creation of a function: needs a path passed to it
  #literally just checking if the files have the same number of ids/seq
  get_read_names <- function(path){
    #from ShortRead: really convinient to just read in fasta
    seqs <- readFastq(path)
    #replaces all the seq ids as characters, then replaces everything in array with empty ""
    sub(as.character(seqs@id), pattern = ".+$", replacement ="")
  }
  #if the expression above is not true, it will return a error message
  #so if the forward and reverse reads are not in matching order
  #comparing the get read names at index 1 and 2 in the paired_file_paths
  #all is comparing the empty arrays and seeing if there is the same number in them
  stopifnot(all(get_read_names(paired_file_paths[1]) == get_read_names(paired_file_paths[2])))
}

#' Take in user's primers and creates the complement, reverse, reverse complement of primers in one tibble
#'
#' @param primer_path a path to the csv file that holds the primer info for this project
#'
#' @return A tibble that contains the forward and reverse primers with all complements of primers
#' @export
#'
#' @examples primer_data <- prepare_primers(primer_path). Part of the larger 'prepare_reads' function. 
prepare_primers <- function(primer_path){
  primer_data_path <- file.path(primer_path)
  primer_data <- read_csv(primer_data_path)
  
  #seperate forward and reverse to make various primers
  forward_primers <- primer_data[, c(1:2)]
  toString(forward_primers)
  reverse_primers <- primer_data[, c(1,3)]
  
  forward_primers$f_compt <- map_chr(forward_primers$forward, function(x) toString(Biostrings::complement(DNAString(x))))
  forward_primers$f_rev <- map_chr(forward_primers$forward, function(x) toString(Biostrings::reverse(DNAString(x))))
  forward_primers$f_rc <- map_chr(forward_primers$forward, function(x) toString(Biostrings::reverseComplement(DNAString(x))))
  
  reverse_primers$r_compt <- map_chr(reverse_primers$reverse, function(x) toString(Biostrings::complement(DNAString(x))))
  reverse_primers$r_rev <- map_chr(reverse_primers$reverse, function(x) toString(Biostrings::reverse(DNAString(x))))
  reverse_primers$r_rc <- map_chr(reverse_primers$reverse, function(x) toString(Biostrings::reverseComplement(DNAString(x))))
  
  #add back together
  primer_data <- forward_primers %>%
    left_join(reverse_primers, by = "primer_name")
  
  return(primer_data)
}

#' A function for calling read_fastq, matching_order_primer_check, and remove_ns functions. This will process and edit the FASTQ and make them ready for the trimming of primers with Cutadapt. Part of a larger 'prepare_reads' function.  
#' @inheritParams dada2::filterAndTrim
#' @param raw_path A path to a directory that contains raw data
#' @param intermediate_path A path to the intermediate folder and directory
#' @param metadata A metadata containing the concatenated metadata and primer data

#'
#' @return
#' @export #add params
#'
#' @examples  Part of the function prepare_reads. fastq_data <- prepare_fastq(raw_path, intermediate_path, maxN = maxN, multithread = multithread)
prepare_fastq <- function(raw_path,intermediate_path, maxN= 0, multithread = FALSE){
  fastq_data <- read_fastq(raw_path)
  matching_order_primer_check(fastq_data)
  fastq_data <- remove_ns(fastq_data, intermediate_path, maxN, multithread=multithread)
  
  return(fastq_data)
}

#Prepare reads. A wrapper function to prepare reads for trimming using Cutadapt. Counts of primers on reads will be output.
#' Main command prepare reads for primer trimming
#' @param directory_path A path to a directory containing reads and metadata/primer files
#' @param primer_path The primer data tibble created in prepare_primers function
#' @param metadata_path A path to a metadata containing the concatenated metadata and primer data
#' @param fastq_path path to a directory containing FASTQ reads
#' @param intermediate_path A path to the intermediate folder and directory
#' @inheritParams prepare_fastq
#'
#' @return A list containing data modified by cutadapt, primer data, FASTQ data, and concatenated metadata and primer data
#' @export
#'
#' @examples returnList<-prepare_reads(directory_path, primer_path, metadata_path, fastq_path, intermediate_path, maxN=0, multithread=TRUE)

prepare_reads <- function(directory_path, primer_path, metadata_path, fastq_path,intermediate_path, maxN=0, multithread=FALSE){
  prep_tables <- file.path(intermediate_path, "Prep_tables.Rdata")
  primer_data <- prepare_primers(primer_path)
  metadata <- prepare_metadata(metadata_path, primer_data)
  fastq_data <- prepare_fastq(raw_path, intermediate_path, maxN = maxN, multithread = multithread)
  pre_primer_hit_data<- pre_primer_hit_data(primer_data, fastq_data,
                                            intermediate_path)
  pre_primer_plot <- primer_hit_plot(pre_primer_hit_data, fastq_data,
                                     intermediate_path, "pre_primer_plot.pdf")
  cutadapt_data <- cutadapt_tibble(fastq_data, metadata, intermediate_path)
  returnList <- list(cutadapt_data=cutadapt_data, primer_data=primer_data, fastq_data=fastq_data, metadata=metadata)
  return(returnList)
}