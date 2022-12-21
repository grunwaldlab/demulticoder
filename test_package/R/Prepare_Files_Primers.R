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
  metadata <- metadata[order(metadata$sample_name),]
  if ("primer_name" %in% colnames(metadata)) { # Check if primer_name is in metadata columns, if true, moves it to the first column
    new_metadata_cols <- c("primer_name", colnames(metadata)[colnames(metadata) != "primer_name"]) 
    metadata = metadata[new_metadata_cols]
    metadata$sample_nameBarcode<-paste0(metadata$sample_name,"_",metadata$primer_name)
    metadata <- metadata %>% relocate(sample_nameBarcode, .after=sample_name)
  }
  else { # Raises error if false
    stop("Please make sure that there is a 'primer_name' column in your metadata table.")
  }
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
                       sample_name = sub(fastq_paths, pattern = "_.+$", replacement = ""),
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
primer_check <- function(fastq_data){
  paired_file_paths <- fastq_data %>%
    filter(sample_name == gdata::first(sample_name)) %>%
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

#' A function for calling read_fastq, primer_check, and remove_ns functions. This will process and edit the FASTQ and make them ready for the trimming of primers with Cutadapt. Part of a larger 'prepare_reads' function.  
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
  primer_check(fastq_data)
  fastq_data <- remove_ns(fastq_data, intermediate_path, maxN, multithread=multithread)
  
  return(fastq_data)
}

#' Get primer counts fo reach sample before primer removal and trimming steps
#'
#' @param primer_data The primer data tibble created in prepare_primers function
#' @param fastq_data A tibble with the fastq file paths, the direction of the sequences, and names of sequences
#' @param intermediate_path A path to the intermediate folder and directory
#'
#' @return A number of reads in which the primer is found
#' @export
#'
#' @examples 
get_pre_primer_hits<- function(primer_data, fastq_data, intermediate_path){
  primer_hit_data <- gather(primer_data, key = "orientation", value = "sequence", forward, f_compt, f_rev,
                            f_rc, reverse, r_compt, r_rev, r_rc)
  
  #from DADA2
  primer_hits <- function(primer, path) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(path)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  
  primer_hit_data_csv_path <- file.path(intermediate_path, "primer_hit_data_pretrim.csv")
  #if a file exists in there, then write the path to it
  if (file.exists(primer_hit_data_csv_path)){
    primer_hit_data <- read_csv(primer_hit_data_csv_path)
    print(primer_hit_data)
    #if the file doesn't exist in primer hit data
    #map applied a function to each part of that specific vector
  } else {
    primer_hit_counts <- future_map(fastq_data$prefiltered_path,
                                    function (a_path) map_dbl(primer_hit_data$sequence, primer_hits, path = a_path))
    #gets or sets the name of an object
    names(primer_hit_counts) <- paste0(fastq_data$file_id)
    #primer hit data will be a tibble with the columns of primer hit data plus the primer hit counts
    #with names pulled from the fastq_data file
    primer_hit_data <- bind_cols(primer_hit_data, as_tibble(primer_hit_counts))
    print(primer_hit_data)
    write_csv(primer_hit_data, primer_hit_data_csv_path)
  }
  return(primer_hit_data)
  
}

#' Get primer counts fo reach sample after primer removal and trimming steps
#'
#' @param primer_data The primer data tibble created in prepare_primers function
#' @param cutadapt_data Intermediate_data folder with trimmed and filtered reads for each sample
#' @param intermediate_path A path to the intermediate folder and directory
#'
#' @return
#' @export
#'
#' @examples
get_post_primer_hits <- function(primer_data, cutadapt_data, intermediate_path){
  post_primer_hit_data <- gather(primer_data, key = "orientation", value = "sequence", forward, f_compt, f_rev,
                                 f_rc, reverse, r_compt, r_rev, r_rc)
  
  #from DADA2
  post_primer_hits <- function(primer, path) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(path)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  
  
  post_primer_hit_data_csv_path <- file.path(intermediate_path, "primer_hit_data_post_trim.csv")
  #if a file exists in there, then write the path to it
  if (file.exists(post_primer_hit_data_csv_path)){
    post_primer_hit_data <- read_csv(post_primer_hit_data_csv_path)
    #if the file doesn't exist in primer hit data
    #map applied a function to each part of that specific vector
  } else {
    post_primer_hit_counts <- future_map(cutadapt_data$filtered_path,
                                         function (a_path) map_dbl(post_primer_hit_data$sequence, post_primer_hits, path = a_path))
    #gets or sets the name of an object
    names(post_primer_hit_counts) <- paste0(cutadapt_data$file_id, "_", cutadapt_data$primer_name)
    #primer hit data will be a tibble with the columns of primer hit data plus the primer hit counts
    #with names pulled from the fastq_data file
    post_primer_hit_data <- bind_cols(post_primer_hit_data, as_tibble(post_primer_hit_counts))
    print(post_primer_hit_data)
    write_csv(post_primer_hit_data, post_primer_hit_data_csv_path)
  }
  return(post_primer_hit_data)
}

#' Make a barplot of primers identified on reads
#'
#' @param primer_hits A number of reads in which the primer is found
#' @param fastq_data A tibble with the fastq file paths, the direction of the sequences, and names of sequences
#' @param plot_name A filename under which a PDF file of the plot will be saved as
#'
#' @return
#' @export
#'
#' @examples
make_primer_hit_plot <- function(primer_hits, fastq_data, intermediate_path, plot_name){
  #This function takes so long to run - wonder if it will be ok with the servers - more efficient way to do this?
  #removing the sequence in the primer_hits tibble
  primer_hits <- primer_hits[-(3)]
  #concatonating the two beginning columns
  primer_hits$primer_type <- paste(primer_hits$primer_name, primer_hits$orientation)
  #subsetting just the concatonated column
  new_primer_hits <- primer_hits[-(1:2)]
  #move name to the first column
  new_primer_hits <- new_primer_hits %>%
    select(primer_type, everything())
  
  #easier to work without the character names
  only_counts <- new_primer_hits[-(1)]
  #add a prefix to the columns so they are easier to mutate
  colnames(only_counts) <- paste("Col", colnames(only_counts), sep="-")
  
  #add them all up
  total_nums <- only_counts %>%
    rowwise() %>%
    mutate(Total= sum(across(starts_with("Col")), na.rm = TRUE))
  
  #add back to new_primer_hits
  new_primer_hits$Total <- paste(total_nums$Total)
  
  #subset just the names and the totals
  needed_cols <- c("primer_type", "Total")
  total_primers <- new_primer_hits[needed_cols]
  
  #convert the total column to numeric
  total_primers <- transform(total_primers, Total = as.numeric(Total))
  #a bar chart
  plot <- ggplot(data=total_primers, aes(x=primer_type, y=Total)) +
    geom_bar(stat="identity", width=1.0, fill= "seagreen3") +
    geom_text(aes(label=Total)) +
    coord_flip()
  print(plot)
  ggsave(plot, filename = plot_name, path = intermediate_path, width = 8, height = 8)
  return(plot)
  
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
  pre_primer_hit_data<- get_pre_primer_hits(primer_data, fastq_data,
                                            intermediate_path)
  pre_primer_plot <- make_primer_hit_plot(pre_primer_hit_data, fastq_data,
                                     intermediate_path, "pre_primer_plot.pdf")
  cutadapt_data <- make_cutadapt_tibble(fastq_data, metadata, intermediate_path)
  returnList <- list(cutadapt_data=cutadapt_data, primer_data=primer_data, fastq_data=fastq_data, metadata=metadata)
  return(returnList)
}