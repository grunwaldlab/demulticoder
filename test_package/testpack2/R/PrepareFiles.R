#' Read in the metadata from user and combine it with the primer data
#'
#' @param metadata_path Path to the metadata csv
#' @param primer_data The primer data tibble created in prepare_primers function
#'
#' @return A metadata containing the concatenated metadata and primer data
#' @export
#'
#' @examples
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

#' Create intermediate path and directory
#'
#' @param working_dir_path
#'
#' @return intermediate_path
#' @export
#'
#' @examples
create_intermediate <- function(working_dir_path){
  intermediate_read_dir <- file.path(working_dir_path,  "intermediate_data")
  if (! dir.exists(intermediate_read_dir)){
    dir.create(intermediate_read_dir)
  }
  intermediate_path <- file.path(intermediate_read_dir)
  return(intermediate_path)
}

#' Takes in the fastq files from the user and creates a tibble with the paths to future used files
#'
#' @param raw_path The path to the fastq files from user
#'
#' @return A tibble with the fastq file paths, the direction of the sequences, and names of sequences
#' @export
#'
#' @examples
read_fastq <- function(raw_path){
  if (BiocGenerics::grepl(fastq_paths, pattern = "_R2" | "_2") == FALSE) {
    stop("Please ensure that the FASTQ files for a sample's read 2 is in the same folder as read 1.")
  }
  
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
#' @param fastq_data
#'
#' @return None
#' @export
#'
#' @examples
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
