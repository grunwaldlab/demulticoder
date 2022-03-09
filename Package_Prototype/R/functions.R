#' Prepare Primers 
#'
#' Take in user primers and creates the complement, reverse, reverse complement of primers in one tibble
#' @param primer_path a path to the csv that holds the primer info for this project.
#' @return A tibble that contains the forward and reverse primers with created complement, etc. 
#' @export
prepare_primers <- function(primer_path){
  primer_data_path <- file.path(primer_path)
  primer_data <- read_csv(primer_data_path)
  
  #seperate forward and reverse to make various primers
  forward_primers <- primer_data[, c(1:2)]
  reverse_primers <- primer_data[, c(1,3)]
  
  
  forward_primers$f_compt <- map_chr(forward_primers$forward, function(x) toString(complement(DNAString(x))))
  forward_primers$f_rev <- map_chr(forward_primers$forward, function(x) toString(reverse(DNAString(x))))
  forward_primers$f_rc <- map_chr(forward_primers$forward, function(x) toString(reverseComplement(DNAString(x))))
  
  reverse_primers$r_compt <- map_chr(reverse_primers$reverse, function(x) toString(complement(DNAString(x))))
  reverse_primers$r_rev <- map_chr(reverse_primers$reverse, function(x) toString(reverse(DNAString(x))))
  reverse_primers$r_rc <- map_chr(reverse_primers$reverse, function(x) toString(reverseComplement(DNAString(x))))
  
  #add back together 
  primer_data <- forward_primers %>%
    left_join(reverse_primers, by = "primer_name")
  
  print(primer_data)
  
  return(primer_data)
  
}

#' Prepare Metadata 
#'
#' Read in the metadata from user and combine it with the primer data
#' @param metadata_path Path to the metadata csv
#' @param primer_data The primer data tibble created in prepare_primers function
#' @return A tibble containing the concatenated metadata and primer data
#' @export
prepare_metadata <- function(metadata_path, primer_data){
  metadata <- read_csv(metadata_path) %>%
    left_join(primer_data, by = c("primer_name"))
  #this could cause issues in other data sets, if the name of the sample is no the first column
  names(metadata)[1] <- "sample_id"
  return(metadata)
}

#' Create Intermediate 
#'
#' Creates a file in the current working directory used for program files later on
#' @return None
#' @export
create_intermediate <- function(){
  if (! dir.exists("intermediate_data")) {
    dir.create("intermediate_data")
  }
}

#' Read Fastq
#'
#' Takes in the fastq files from the user and creates a tibble with the paths to future used files
#' @param fastq_path The path to the fastq files from user
#' @return A tibble with the fastq file paths, the direction of the sequences, and names of sequences
#' @export
read_fastq <- function(fastq_path){
  fastq_paths <- list.files(fastq_path, pattern = "\\.fastq")
  #constructing a tibble - special data frame with improved behaviors. 
  fastq_data <- tibble(file_id = sub(fastq_paths, pattern = "\\.fastq\\.gz$", replacement = ""),
                       sample_id = sub(fastq_paths, pattern = "_.+$", replacement = ""),
                       #direction: the grep1 is looking to match the pattern, relates a TRUE or FALSE with it (aka 1 or 0) and adds 1 to find the 
                       #correct index to put Reverse or Forward
                       direction = c( "Reverse", "Forward")[BiocGenerics::grepl(fastq_paths, pattern = "_1") + 1],
                       raw_path = file.path(fastq_path, fastq_paths))
  
  return(fastq_data)
}

#' Mathcing Order Primer Check
#'
#' Validation for the primers - making sure the forward and reverse reads are in matching order
#' @param fastq_data The tibble that the read_fastq function spits out
#' @return None
#' @export
matching_order_primer_check <- function(fastq_data){
  paired_file_paths <- fastq_data %>%
    filter(sample_id == first(sample_id)) %>%
    #pull is just like $ - you are accessing a variable inside the tibble 
    pull(raw_path)
  
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

#' Remove N's
#'
#' Removes any sequences that have Ns in them
#' @param fastq_data The tibble that the read_fastq function spits out
#' @param intermediate_path The path to the intermediate data file 
#' @return A tibble with the sequences with the N's removed 
#' @export
remove_ns <- function(fastq_data, intermediate_path){
  prefiltered_read_dir <- file.path(intermediate_path, "prefiltered_sequences")
  fastq_data$prefiltered_path <- file.path(prefiltered_read_dir, base::basename(fastq_data$raw_path))
  #if the files do not exist in the prefiltered path path (which clearly they don't)
  #raw path is the path to the actual fastq files in mock community 
  #fwd takes from the raw path data 
  #filt is the path to the output filtered files that we created 
  #the sequences it chose to put into prefiltered path had to have no Ns in them 
  if(! all(file.exists(fastq_data$prefiltered_path))){
    filterAndTrim(fwd = fastq_data[fastq_data$direction == "Forward", ][["raw_path"]],
                  filt = fastq_data[fastq_data$direction == "Forward", ][["prefiltered_path"]],
                  rev = fastq_data[fastq_data$direction == "Reverse", ][["raw_path"]],
                  filt.rev = fastq_data[fastq_data$direction == "Reverse", ][["prefiltered_path"]],
                  maxN = 0,
                  multithread = TRUE)
  }
  
  return(fastq_data)
  
}

#' Prepare Fastq
#'
#' Does all of the processing/editing of the fastqs in order to make ready for cutadapt
#' @param fastq_path The path to the fastq files
#' @param intermediate_path The path to the intermediate data file 
#' @return A data frame with the fastq information 
#' @export
prepare_fastq <- function(fastq_path,intermediate_path){
  fastq_data <- read_fastq(fastq_path)
  matching_order_primer_check(fastq_data)
  #this takes a hot sec
  fastq_data_filtered <- remove_ns(fastq_data, intermediate_path)
  
  #crate a place to put the results of both reads that were trimmed successfully and those that were not 
  trimmed_read_dir <- file.path(intermediate_path, "trimmed_sequences")
  if (! dir.exists(trimmed_read_dir)){
    dir.create(trimmed_read_dir)
  }
  
  fastq_data_filtered$trimmed_path <- file.path(trimmed_read_dir, paste0(fastq_data_filtered$file_id, ".fastq.gz"))
  
  untrimmed_read_dir <- file.path(intermediate_path, "untrimmed_sequences")
  if (! dir.exists(untrimmed_read_dir)) {
    dir.create(untrimmed_read_dir)
  }
  fastq_data_filtered$untrimmed_path <- file.path(untrimmed_read_dir, paste0(fastq_data_filtered$file_id, ".fastq.gz"))
  
  return(fastq_data_filtered)
}

#' Pre Primer Hit Data
#'
#' To give a visualization (in a table) showing the different primer hits in each of the sequences for each sample
#' @param primer_data The primer data tibble created in prepare_primers function
#' @param fastq_data The tibble that the read_fastq function spits out
#' @param intermediate_path The path to the intermediate data file 
#' @return A data frame with the primer hits in each sample
#' @export
pre_primer_hit_data <- function(primer_data, fastq_data, intermediate_path){
  primer_hit_data <- gather(primer_data, key = "orientation", value = "sequence", forward, f_compt, f_rev, 
                            f_rc, reverse, r_compt, r_rev, r_rc)
  
  #from DADA2
  primer_hits <- function(primer, path) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(path)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  
  
  primer_hit_data_csv_path <- file.path(intermediate_path, "primer_hit_data_test.csv")
  #if a file exists in there, then write the path to it 
  if (file.exists(primer_hit_data_csv_path)){
    primer_hit_data <- read_csv(primer_hit_data_csv_path)
    #if the file doesn't exist in primer hit data 
    #map applied a function to each part of that specific vector 
  } else {
    primer_hit_counts <- future_map(fastq_data$prefiltered_path, 
                                    function (a_path) map_dbl(primer_hit_data$sequence, primer_hits, path = a_path))
    #gets or sets the name of an object 
    names(primer_hit_counts) <- fastq_data$file_id
    #primer hit data will be a tibble with the columns of primer hit data plus the primer hit counts 
    #with names pulled from the fastq_data file
    primer_hit_data <- bind_cols(primer_hit_data, as_tibble(primer_hit_counts))
    write_csv(primer_hit_data, primer_hit_data_csv_path)
  }
  return(primer_hit_data)
}

#' Primer Hit Plot
#'
#' To create a plot showing the primer hits for a sanity check
#' @param primer_hits A tibble that shows the primer hits, produced by pre_primer_hit_data
#' @param fastq_data The tibble that the read_fastq function spits out
#' @param metadata The metadata table that prepare_metadata provides
#' @return A plot that shows the primer hits 
#' @export
primer_hit_plot <- function(primer_hits, fastq_data, metadata){
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
  
  return(plot)
}


#################################################################################################################
## Name: change_sample_ids
## Purpose: THIS IS NOT A ROBUST FUNCTION - NEEDED IT QUICKLY
## Inputs:
## Returns:
##################################################################################################################
change_sample_ids <- function(fastq_data){
  #change the sample id's in the fastq file (ONLY for Zach's data)
  library_data <- c("A1", "A3")
  sample_data <- c("SRR13658662", "SRR13658661")
  overall <- data.frame(library_data, sample_data)
  
  fastq_data2 <- fastq_data %>%
    right_join(overall, by = c("sample_id" = "sample_data")) %>%
    select(-sample_id) %>%
    rename(sample_id = library_data) %>%
    relocate(sample_id, .before = direction)
  
  return(fastq_data2)
}


#' Cutadapt Tibble
#'
#' Creates a tibble that can be fed into cutadapt by combining metadata, primer, and fastq tibbles
#' @param new_fastq_data The fastq data that came out of prepare fastq
#' @param metadata The metadata table that prepare_metadata provides
#' @return A tibble for cutadapt
#' @export
cutadapt_tibble <- function(fastq_data, metadata){
  cutadapt_data <- metadata %>%
    select(sample_id, forward, f_compt, f_rev, f_rc, reverse, r_compt, r_rev, r_rc)
  
  cutadapt_data <- new_fastq_data %>%
    #choosing only forward from the fastq data file 
    filter(direction == "Forward") %>%
    select(sample_id, prefiltered_path, trimmed_path, untrimmed_path) %>%
    right_join(cutadapt_data, by = "sample_id") %>%
    rename(forward_input_path = prefiltered_path,
           forward_output_path = trimmed_path, 
           forward_untrimmed_path = untrimmed_path) 
  #chossing only reverse from the fastq data file
  cutadapt_data <- new_fastq_data %>%
    filter(direction == "Reverse") %>%
    select(sample_id, prefiltered_path, trimmed_path, untrimmed_path) %>%
    right_join(cutadapt_data, by = "sample_id") %>%
    rename(reverse_input_path = prefiltered_path,
           reverse_output_path = trimmed_path,
           reverse_untrimmed_path = untrimmed_path)
  
  return(cutadapt_data)
}

#' Post Primer Hit Data
#'
#' Makes a primer hit tibble for fastq sequences AFTER cutadapt has been run
#' @param primer_hit_data the hit data from pre_primer_hit data 
#' @param primer_data The tibble from prepare_primers
#' @param fastq_data The tibble that the read_fastq function spits out
#' @param intermediate_path The path to the intermediate data file 
#' @return A tibble with all the primer hits 
#' @export
post_primer_hit_data <- function(primer_hit_data, primer_data, fastq_data, intermediate_path){
  #DADA2 function...
  primer_hits <- function(primer, path) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(path)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  
  cutadapt_verify_path <- file.path(intermediate_path, "cutadapt_verify_data.csv")
  #if a file exists in there, then write the path to it 
  if (file.exists(cutadapt_verify_path)){
    cutadapt_verify_data <- read_csv(cutadapt_verify_path)
    #if the file doesn't exist in primer hit data 
    #map applied a function to each part of that specific vector 
  } else {
    primer_hit_counts <- future_map(fastq_data$trimmed_path, 
                                    function (a_path) map_dbl(primer_hit_data$sequence, primer_hits, path = a_path))
    #gets or sets the name of an object 
    names(primer_hit_counts) <- fastq_data$file_id
    #primer hit data will be a tibble with the columns of primer hit data plus the primer hit counts 
    #with names pulled from the fastq_data file
    cutadapt_verify_data <- primer_hit_data
    cutadapt_verify_data[names(primer_hit_counts)] <- primer_hit_counts
    write_csv(cutadapt_verify_data, cutadapt_verify_path)
    
  }
  return(cutadapt_verify_data)
}

#' Filter and Trim
#'
#' To filter reads after cutadapt with various choices - i.e. with Ns, below a length of 50, or greater than expected error limit
#' @param intermediate_path The path to the intermediate data file 
#' @param fastq_data The tibble that the read_fastq function spits out
#' @return None
#' @export
filter_and_trim <- function(intermediate_path, fastq_data){
  #make a place to put the results 
  expected_error_filter_limit <- 5 #default setting
  truncation_qual_limit <- 5 #default setting 
  filtered_reads_dir <- file.path(intermediate_path, "filtered_sequences")
  fastq_data$filtered_path <- file.path(filtered_reads_dir, paste0(fastq_data$file_id, ".fastq.gz"))
  
  if(! all(file.exists(fastq_data$filtered_path))){
    filter_results <- filterAndTrim(fwd = fastq_data$trimmed_path[fastq_data$direction == "Forward"], 
                                    filt = fastq_data$filtered_path[fastq_data$direction == "Forward"],
                                    rev =  fastq_data$trimmed_path[fastq_data$direction == "Reverse"], 
                                    filt.rev = fastq_data$filtered_path[fastq_data$direction == "Reverse"], 
                                    maxN = 0, 
                                    maxEE = c(expected_error_filter_limit, expected_error_filter_limit), 
                                    truncQ = truncation_qual_limit, 
                                    minLen = 50, 
                                    rm.phix = TRUE, 
                                    compress = TRUE,
                                    multithread = TRUE)
    filter_results <- as_tibble(filter_results)
    print(colMeans(filter_results))
    
  }
}

#' Cutadapt Run
#'
#' To check for the user having cutadapt - if yes: run the cutadapt command
#' @param cutadapt_path The path to the acutal cutadapt program
#' @param cutadapt_data The tibble from cutadapt tibble function
#' @return None
#' @export
cutadapt_run <- function(cutadapt_path, cutadapt_data){
  cutadapt <- cutadapt_path
  tryCatch(system2("cutadapt", args = "--version"),
           warning=function(w){
             stop("cutadapt cannot be found on PATH. Check if it is installed?")
           })
  #R1_flags: -g is the 5' adaptor (the forward primer) and -a is the reverse
  #make a flag for both forward and reverse reads 
  #-o is the forward_output_path (where we want cutadapt to store the trimming)
  #-n is by defalt 1 (by default only one adapter sequence is removed from each read, but we 
  #want 2, so have to specify it)
  #-p is the reverse_output_path (where we want cutadapt to store the trimming)
  #--minimum length = the min length reads have to be to be kept
  
  
  #we are jsut building a column in cutadapt data that has the command line arguements in it
  cutadapt_data <- cutadapt_data %>%
    mutate(command_args = paste(
      "-g", forward,
      "-a", r_rc,
      "-G", reverse,
      "-A", f_rc,
      "-n", 2,
      "-o", forward_output_path,
      "-p", reverse_output_path,
      "--minimum-length", 50,
      "--untrimmed-output", forward_untrimmed_path,
      "--untrimmed-paired-output", reverse_untrimmed_path,
      "--quiet",
      forward_input_path,
      reverse_input_path
    ))
  
  #check if the commands were correctly entered 
  print(cutadapt_data)
  
  #now actually run the cutadapt 
  #should check if all the files for cutadapt even exist, but thats a future problem TODO
  
  cutadapt_output <- future_map(cutadapt_data$command_args, ~system2("cutadapt", args = .x))
  
}