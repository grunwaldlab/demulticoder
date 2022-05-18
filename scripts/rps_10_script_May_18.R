library(dada2)
library(ShortRead)
library(Biostrings)
library(dplyr)
library(purrr)
library(furrr)
library(tidyr)
library(readr)
library(ggplot2)
library(gridExtra)
library(metacoder)
library(stringr)
library(sessioninfo)

seed <- 1
set.seed(seed)
##make sure people have the proper versions installed-way to make this streamlined
##TODO
##R studio format for easier viewing of pipeline

#****************************************************************FUNCTIONS TO ANALYZE METABARCODING DATA (IN PROGRESS)********************************
####TODO: 
# Further generalize functions
# Make a main function(s) (the ones the users would actually get to interact with)
# Markdown file for documentation
# Error catching in the functions 
#****************************************************************FUNCTION DECLARATIONS*****************************************************************

#################################################################################################################
## Name: prepare_primers
## Purpose: take in user primers and creates the complement, reverse, reverse complement of primers in one tibble
## Inputs: a path to the csv that holds the primers needed for project. Particular format needed
## Returns: a tibble that contains the forward and reverse primers with created complement, etc. 
################################################################################################################## 
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
  
  print(primer_data)
  
  return(primer_data)
  
}


#################################################################################################################
## Name: prepare_metadata
## Purpose: to read in the metadata from user and combine it with the primer data
## Inputs: the metadata path, primer data tibble created in prepare_primers()
## Returns: a data frame containing the metadata and primer data now combined
##################################################################################################################
prepare_metadata <- function(metadata_path, primer_data){
  metadata <- read_csv(metadata_path) %>%
    left_join(primer_data, by = c("primer_name"))

  # (Hung) error handling, if there is a "primer_name" column, move it to the first column. If not, prints error
  if ("primer_name" %in% colnames(metadata)) {
    new_metadata_cols <- c("primer_name", colnames(metadata)[colnames(metadata) != "primer_name"])
    metadata = metadata[new_metadata_cols]
  }
  else {
    stop("Please make sure that there is a 'primer_name' column in your metadata table.")
  }
  #this could cause issues in other data sets, if the name of the sample is no the first column
  # names(metadata)[1] <- "sample_id"
  ## (Hung) The line above moves the 'sample_id' column to the first column in the data frame. I don't want any complication with the control flow.
  return(metadata)
}


#################################################################################################################
## Name: create_intermediate
## Purpose: to create a file in the current working directory used for program files later on
## Inputs: working directory path
## Returns: Makes an intermediate data folder
##################################################################################################################
create_intermediate <- function(working_dir_path){
  intermediate_read_dir <- file.path(working_dir_path,  "intermediate_data")
  if (! dir.exists(intermediate_read_dir)){
    dir.create(intermediate_read_dir)
  }
}


#################################################################################################################
## Name: read_fastq 
## Purpose: take in the fastq files from the user and creating a tibble that organizes them
## Inputs: fastq file path 
## Returns: a tibble with the fastq file paths, the direction of the sequences, and names of sequences
##################################################################################################################
read_fastq <- function(raw_path){
  fastq_paths <- list.files(raw_path, pattern = "\\.fastq")
  #constructing a tibble - special data frame with improved behaviors. 
  fastq_data <- tibble(file_id = sub(fastq_paths, pattern = "\\.fastq\\.gz$", replacement = ""),
                       sample_id = sub(fastq_paths, pattern = "_.+$", replacement = ""),
                       #direction: the grep1 is looking to match the pattern, relates a TRUE or FALSE with it (aka 1 or 0) and adds 1 to find the 
                       #correct index to put Reverse or Forward
                       direction = c("Reverse", "Forward")[BiocGenerics::grepl(fastq_paths, pattern = "_R1" | "_1") + 1] | # (Hung) | is the OR operator, returns true if either of these exists
                       raw_data_path = file.path(raw_path, fastq_paths))
  
  # (Hung) Check if the fastq_paths have files with _R2 or _2 pattern in file name. If false, raises error.
  if (BiocGenerics::grepl(fastq_paths, pattern = "_R2" | "_2") == FALSE) {
    stop("Please ensure that the FASTQ files for a sample's read 2 is in the same folder as read 1.")
  }
  
  return(fastq_data)
}



#################################################################################################################
## Name: matching_order_primer_check
## Purpose: validation for the primers - making sure the forward and reverse reads are in matching order
## Inputs: the fastq data tibble
## Returns: none
##################################################################################################################
matching_order_primer_check <- function(fastq_data){
  paired_file_paths <- fastq_data %>%
    filter(sample_id == first(sample_id)) %>%
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

#################################################################################################################
#################################################################################################################
## Name: remove_ns
## Purpose: removes any sequences that have Ns in them
## Inputs: the fastq data tibble, path to the intermediate data 
## Returns: TODO return the data frame with the fastq_data
##################################################################################################################
remove_ns <- function(fastq_data, intermediate_path, metadata){
  #fastq_data <- metadata %>%
  #left_join(fastq_data, by = c("sample_id" = "sample_id"))
  prefiltered_read_dir <- file.path(intermediate_path, "prefiltered_sequences")
  fastq_data$prefiltered_path <- file.path(prefiltered_read_dir, base::basename(fastq_data$raw_data_path))
  #if the files do not exist in the prefiltered path path (which clearly they don't)
  #raw path is the path to the actual fastq files in mock community 
  #fwd takes from the raw path data 
  #filt is the path to the output filtered files that we created 
  #the sequences it chose to put into prefiltered path had to have no Ns in them 
  if(!all(file.exists(fastq_data$prefiltered_path))){
    filterAndTrim(fwd = fastq_data[fastq_data$direction == "Forward", ][["raw_data_path"]],
                  filt = fastq_data[fastq_data$direction == "Forward", ][["prefiltered_path"]],
                  rev = fastq_data[fastq_data$direction == "Reverse", ][["raw_data_path"]],
                  filt.rev = fastq_data[fastq_data$direction == "Reverse", ][["prefiltered_path"]],
                  maxN = 0,
                  multithread = TRUE)
  }
  
  return(fastq_data)
  
}
#################################################################################################################
## Name: prepare_fastq (don't run above function individually)
## Purpose: does all of the processing/editing of the fastqs in order to make ready for cutadapt
## Inputs: path to the fastq files, the intermediate data path
## Outputs: a data frame with the fastq information 
##################################################################################################################
prepare_fastq <- function(raw_path,intermediate_path, metadata){
  fastq_data <- read_fastq(raw_path)
  matching_order_primer_check(fastq_data)
  fastq_data <- remove_ns(fastq_data, intermediate_path, metadata)
  
  
  return(fastq_data)
}


#removed tibble merging steps, and will do this later during Cutadapt prep steps since it creates complications when I run the primer hit functions

#################################################################################################################
## Name: pre_primer_hit_data
## Purpose: To give a visualization (in a table) showing the different primer hits in each of the sequences for 
## each sample
## Inputs: primer data tibble (from prepare_primers()), fastq data tibble (from prepare_fastq()), and the intermediate
## data path
## Returns: A tibble of the primer hits
##################################################################################################################
pre_primer_hit_data <- function(primer_data, fastq_data, intermediate_path){
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
    write_csv(primer_hit_data, primer_hit_data_csv_path)
  }
  return(primer_hit_data)
}

#################################################################################################################
## Name: primer_hit_plot
## Purpose: To create a plot showing the primer hits for a sanity check
## Inputs: the primer hit tibble, fastq_data, metadata
## Returns: a plot that shows the primer hits 
##################################################################################################################
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
  ggsave(plot, filename = 'primer_hits.pdf', path = 'intermediate_data', width = 8, height = 8)
  
}


#Adjust file name based on function info

#ggsave(merge_plot, filename = 'read_merging.png', path = 'results', width = 8, height = 8)
#TODO save output plot into new output folder with results

#################################################################################################################
## Name: cutadapt_tibble
## Purpose: create a tibble that can be fed into cutadapt 
## Inputs: a fastq_data tibble created by prepare_fastq
## Returns: a cutadapt tibble
##################################################################################################################
cutadapt_tibble <- function(fastq_data, metadata, intermediate_path){ #new_fastq_data needed, why not just fastq_data
  cutadapt_data <- metadata %>%
    left_join(fastq_data, by = c("sample_id" = "sample_id"))
  #remove extra columns with NA
  cutadapt_data <- na.omit(cutadapt_data)
  trimmed_read_dir <- file.path(intermediate_path, "trimmed_sequences")
  if (! dir.exists(trimmed_read_dir)){
    dir.create(trimmed_read_dir)
  }
  
  cutadapt_data$trimmed_path <- file.path(trimmed_read_dir, paste0(cutadapt_data$file_id, "_", cutadapt_data$primer_name, ".fastq.gz"))
  
  #fastq_data_filtered$trimmed_path <- file.path(trimmed_read_dir, paste0(fastq_data_filtered$file_id, ".fastq.gz"))
  untrimmed_read_dir <- file.path(intermediate_path, "untrimmed_sequences")
  if (! dir.exists(untrimmed_read_dir)) {
    dir.create(untrimmed_read_dir)
  }
  cutadapt_data$untrimmed_path <- file.path(untrimmed_read_dir,  paste0(cutadapt_data$file_id, "_", cutadapt_data$primer_name, ".fastq.gz"))
  
  filtered_read_dir <- file.path(intermediate_path, "filtered_sequences")
  if (! dir.exists(filtered_read_dir)) {
    dir.create(filtered_read_dir)
  }
  cutadapt_data$filtered_path <- file.path(filtered_read_dir,  paste0(cutadapt_data$file_id, "_", cutadapt_data$primer_name, ".fastq.gz"))
  
  
  return(cutadapt_data)
}

#################################################################################################################
## Name: cutadapt_run_make it more like Sam's function
## Purpose: to check for the user having cutadapt - if yes: run the cutadapt command
## Inputs: cutadapt path (YOU NEED CUTADAPT PROGRAM), cutadapt data tibble
## Returns: none
##################################################################################################################
#Not exactly like Sam's function. Needs additional simplification and generalization of parameters
cutadapt_run <- function(cutadapt_path, cutadapt_data){
  cutadapt <- cutadapt_path
  tryCatch(system2(cutadapt, args = "--version"),
           warning=function(w){
             stop("cutadapt cannot be found on PATH. Check if it is installed?")
           })
  #simplify
  #may need to demultiplex first-meed to assess this
  cutadapt <- path.expand(cutadapt_path)
  R1_flags = unique(paste("-g",cutadapt_data$forward, "-a", cutadapt_data$r_rc))
  R2_flags = unique(paste("-G",cutadapt_data$reverse, "-A", cutadapt_data$f_rc))
  fwd_trim = cutadapt_data[cutadapt_data$direction == "Forward", ][["trimmed_path"]]
  rev_trim = cutadapt_data[cutadapt_data$direction == "Reverse", ][["trimmed_path"]]
  fwd_untrim = cutadapt_data[cutadapt_data$direction == "Forward", ][["untrimmed_path"]]
  rev_untrim = cutadapt_data[cutadapt_data$direction == "Reverse", ][["untrimmed_path"]]
  fwd_prefilt = cutadapt_data[cutadapt_data$direction == "Forward", ][["prefiltered_path"]]
  rev_prefilt = cutadapt_data[cutadapt_data$direction == "Reverse", ][["prefiltered_path"]]
  #consider parameters
  command_args=paste(R1_flags, R2_flags, "-n", 2,  "-o",  fwd_trim, "-p", rev_trim, 
                      "--minimum-length", 50,"--maximum-length", 290,"--untrimmed-output", fwd_untrim, 
                      "--untrimmed-paired-output", rev_untrim,"--quiet", fwd_prefilt, rev_prefilt)
  if (! all(file.exists(c(cutadapt_data$trimmed_path)))) {
    cutadapt_output <- future_map(command_args, ~system2(cutadapt, args = .x))
  }
  
}
#################################################################################################################
## Name: plot_qc
## Purpose: Function is from DADA2
## Inputs: sample you are interested in plotting (use quotes), cutadapt tibble
## Returns: Plots of quality score at each base position
##################################################################################################################
plot_qc<-function(cutadapt_data, pick_sample){
  #just retrieve all plots for first sample
  sample_info = cutadapt_data$trimmed_path[cutadapt_data$sample_id == pick_sample]
  quality_plots<-plotQualityProfile(sample_info)
  return(quality_plots)
}

#TODO save output plot into new output folder
#ggsave(merge_plot, filename = 'read_merging.png', path = 'results', width = 8, height = 8)

#################################################################################################################
## Name: filter_and_trim
## Purpose: to filter out reads with Ns, below a length of 50, or greater than expected error limis. References filterAndTrim function from DADA2.
## Inputs: intermediate file path, cutadapt tibble
## Returns: none
##################################################################################################################
#Not exactly like Sam's function. Needs additional simplification and generalization of parameters
filter_and_trim <- function(intermediate_path, cutadapt_data){
  #make a place to put the results 
  ## (Hung) Inserted a few variables that can be input by the end user to adjust the parameters of the filterAndTrim function.
  print("Please enter the maximum expected errors (Default: Inf). Reads with more expected errors than this number will be discarded.")
  # expected_error_filter_limit <- 5 #default setting
  expected_error_filter_limit <- readline()
  print("\n")
  print("Please enter a truncation quality limit (Default: 2). Reads at the first instance of a quality score less than or equal to this number will be truncated.")
  # truncation_qual_limit <- 2 #default setting 
  truncation_qual_limit <- as.integer(readline())
  print("\n")
  print("Please enter your desired max_N number (Default: 0). Sequences with more Ns than max_N will be discarded.")
  max_N <- as.integer(readline())
  print("\n")
  print("Please enter your desired min_Len number (Default: 20). Reads containing a quality score of less than min_Len will be discarded.")
  min_Len <- as.integer(readline())
  filtered_read_dir <- file.path(intermediate_path, "filtered_sequences")
  cutadapt_data$filtered_path <- file.path(filtered_read_dir, paste0(cutadapt_data$file_id, "_", cutadapt_data$primer_name, ".fastq.gz"))
  if(! all(file.exists(cutadapt_data$filtered_path))){
    filter_results <- filterAndTrim(fwd = cutadapt_data$trimmed_path[cutadapt_data$direction == "Forward"], 
                                    filt = cutadapt_data$filtered_path[cutadapt_data$direction == "Forward"],
                                    rev =  cutadapt_data$trimmed_path[cutadapt_data$direction == "Reverse"], 
                                    filt.rev = cutadapt_data$filtered_path[cutadapt_data$direction == "Reverse"], 
                                    maxN = max_N, 
                                    maxEE = c(expected_error_filter_limit, expected_error_filter_limit), 
                                    truncQ = truncation_qual_limit, 
                                    minLen = min_Len, 
                                    rm.phix = TRUE, 
                                    compress = TRUE,
                                    matchIDs = TRUE, 
                                    multithread = TRUE)
    filter_results <- as_tibble(filter_results)
    print(colMeans(filter_results))
    
  }
  return(cutadapt_data)
}

#################################################################################################################
## Name: post_primer_hit_data
## Purpose: makes a primer hit tibble for fastq sequences AFTER cutadapt has been run
## Inputs: primer_data information, the fastq_data tibble, the intermediate path
## Returns: a tibble with all the primer hits 
##################################################################################################################
post_primer_hit_data <- function(primer_data, cutadapt_data, intermediate_path){
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
    write_csv(post_primer_hit_data, post_primer_hit_data_csv_path)
  }
  return(post_primer_hit_data)
}

#################################################################################################################
## Name: post_trim_qc
## Purpose: Function is from DADA2
## Inputs: Sample you are interested in plotting (use quotes), cutadapt tibble
## Returns: Plots of quality score at each base position after trimming
##################################################################################################################
post_trim_qc<-function(cutadapt_data, pick_sample){
  #just retrieve all plots for first sample
  sample_info = cutadapt_data$filtered_path[cutadapt_data$sample_id == pick_sample]
  quality_plots2<-plotQualityProfile(sample_info)
  return(quality_plots2)
}

#TODO save output plot into new output folder
#ggsave(merge_plot, filename = 'read_merging.png', path = 'results', width = 8, height = 8)
#################################################################################################################
## Database functions-Feb 2022
## Purpose: Make usable databases that people can reference without having to do significant reformatting during the analysis phase 
## Inputs: Raw databases from http://oomycetedb.cgrb.oregonstate.edu/search.html and https://unite.ut.ee/repository.php
## Returns: Reformatted database and a .csv of genus counts for each database. 
## Misc notes: Still in progress. Adapted from Zach Foster's pipeline scripts. Each database has slightly different header formats, and it is a challenge to generalize these functions. 

##example headers
##Rps10 >name=Aphanomyces_invadans|strain=NJM9701|ncbi_acc=KX405005|ncbi_taxid=157072|oodb_id=13|taxonomy=cellular_organisms;Eukaryota;Stramenopiles;Oomycetes;Saprolegniales;Saprolegniaceae;Aphanomyces_invadans
##ITS Unite header >Cryptosphaeria_sp|KU761147|SH1572218.08FU|reps|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Diaporthales;f__Valsaceae;g__Cryptosphaeria;s__Cryptosphaeria_sp

################################################################################################################## 
#Download and place your desired database into raw_data directory

create_ref_database <- function(intermediate_path){
  formatted_ref_dir <- file.path(intermediate_path, "reference_databases")
  if (! dir.exists(formatted_ref_dir)) {
    dir.create(formatted_ref_dir)
  }
}

format_database <-function(raw_data_path, database_name){
  rps10_db <- read_fasta(file.path(raw_data_path, database_name))
  rps10_data <- str_match(names(rps10_db), pattern = "name=(.+)\\|strain=(.+)\\|ncbi_acc=(.+)\\|ncbi_taxid=(.+)\\|oodb_id=(.+)\\|taxonomy=(.+)$")
  colnames(rps10_data) <- c("header", "name", "strain", "ncbi_acc", "ncbi_taxid", "oodb_id", "taxonomy")
  rps10_data <- as_tibble(rps10_data)
  rps10_data$taxonomy <- gsub(rps10_data$taxonomy, pattern = 'cellular_organisms;', replacement = '', fixed = TRUE)
  rps10_data$taxonomy <- gsub(rps10_data$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  rps10_data$taxonomy <- gsub(rps10_data$taxonomy, pattern = 'Eukaryota', replacement = 'Eukaryota;Heterokontophyta', fixed = TRUE)
  binomial <- map_chr(str_split(rps10_data$taxonomy, pattern = ';'), `[`, 7)
  genus <- map_chr(str_split(binomial, pattern = '_'), `[`, 1)
  unique(genus)
  rps10_data$taxonomy <- map_chr(seq_along(rps10_data$taxonomy), function(index) {
    sub(rps10_data$taxonomy[index], pattern = binomial[index], replacement = paste0(genus[index], ';', binomial[index]))
  })
  rps10_data$taxonomy <- paste0(rps10_data$taxonomy, ';', 'oodb_', seq_along(rps10_data$taxonomy))
  rps10_data$taxonomy <- paste0(rps10_data$taxonomy, ';')
  rps10_data$taxonomy <- trimws(rps10_data$taxonomy)
  rps10_db <- trimws(rps10_db)
  #optional-need to decide on finalized database format
  stopifnot(all(str_count(rps10_data$taxonomy, pattern = ";") == 9))
  genus_count <- table(map_chr(strsplit(rps10_data$name, split = '_'), `[`, 1))
  count_table <- as.data.frame(genus_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Genus', 'Number of sequences')
  formatted_ref_dir <- file.path(intermediate_path, "reference_databases")
  write_csv(count_table, file = file.path(formatted_ref_dir, "rps10_genus_count_table.csv"))
  rps10_ref_path <- file.path(formatted_ref_dir, "rps10_reference_db.fa")
  paste0(">", rps10_data$taxonomy, "\n", rps10_db) %>%
    write_lines(file = rps10_ref_path)
  return(rps10_data)
}

## ITS database
format_database2 <-function(raw_data_path, database_name){
  its_db <- read_fasta(file.path(raw_data_path, database_name))
  its_data <- str_match(names(its_db), pattern = "(.+)\\|(.+)\\|(.+)\\|(.+)\\|(.+)$")
  colnames(its_data) <- c("header", "name", "ncbi_acc", "unite_db", "db", "taxonomy")
  its_data <- as_tibble(its_data)
  its_data$taxonomy <- gsub(its_data$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  its_data$taxonomy <- paste0('Eukaryota;', its_data$taxonomy)
  its_data$taxonomy <- gsub(its_data$taxonomy, pattern = 'Stramenopila;Oomycota', replacement = 'Heterokontophyta;Stramenopiles', fixed = TRUE)
  its_data$taxonomy <- paste0(its_data$taxonomy, ';', 'unite_', seq_along(its_data$taxonomy))
  its_data$taxonomy <- gsub(its_data$taxonomy, pattern = "[a-z]__", replacement = '')
  its_data$taxonomy <- paste0(its_data$taxonomy, ';')
  its_data$taxonomy <- trimws(its_data$taxonomy)
  #Fix after checking out later analysis
  stopifnot(all(str_count(its_data$taxonomy, pattern = ";") == 9))
  genus_count <- table(map_chr(strsplit(its_data$name, split = '_'), `[`, 1))
  count_table <- as.data.frame(genus_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Genus', 'Number of sequences')
  formatted_ref_dir <- file.path(intermediate_path, "reference_databases")
  write_csv(count_table, file = file.path(formatted_ref_dir, "its_genus_count_table.csv"))
  its_ref_path <- file.path(formatted_ref_dir , "its_reference_db.fa")
  paste0(">", its_data$taxonomy, "\n", its_db) %>%
    write_lines(file = its_ref_path)
  return(its_data)
  
}

#create new function for 16S and other potential databases? 

#################################################################################################################
## Infer ASVs-March 2022 STILL INCOMPLETE STARTING AT THIS STEP
## Purpose: Based on existing data/database, make it easier for people to use DADA2 functions to compile ASV abundance tables and taxonomy tables. 
## Inputs: Trimmed reads, Cutadapt tibble, and importantly, formatted databases. (See note above about where to retrieve database fasta files) 
## You need to download the databases first.
## Returns: ASV matrix and taxonomy tables
## Updates: Names of functions could be more clear
################################################################################################################## 
#These are the parameters that Zach used. TODO-Still need to test parameters and modify based on best practices
#min_read_merging_overlap <- 15 # Default is 12
#max_read_merging_mismatch <- 2 # Default is 0
remove_chimeras <- TRUE # Default is TRUE
min_asv_length <- 50
its_clustering_threshold <- 0.99
rps10_clustering_threshold <- 0.96

#################################################################################################################
## Name: Get_fastq_paths; infer_asvs
## Purpose: Core functions are from DADA2, and wrapper function adapted from Zach Foster. Learn error rates and infer ASVs.
## Inputs: Cutadapt data table and path to filtered sequences
## Returns: Merged reads after inference of asvs
##################################################################################################################
#Core DADA2 functions. Pipeline workflow adapted from Zach Foster. DADA2 has some additional parameters as well
#Learn error rates of forward filtered reads
#Select whether you multi-thread or not with TRUE, FALSE

get_fastq_paths <- function(my_direction, my_primer_pair_id) {
  cutadapt_data %>%
    filter(direction == my_direction, primer_name == my_primer_pair_id, file.exists(filtered_path)) %>%
    pull(filtered_path)
}

infer_asvs <-function(my_primer_pair_id, my_direction){
  fastq_paths <- get_fastq_paths(my_direction, my_primer_pair_id)
  error_profile <- learnErrors(fastq_paths, multithread = TRUE)
  
  #if (error_profile) {
  cat(paste0('Error rate plot for the ', my_direction, ' read of primer pair ', my_primer_pair_id, ' \n'))
  plot_errors<-plotErrors(error_profile, nominalQ = TRUE)
  asv_data <- dada(fastq_paths, err = error_profile, multithread = TRUE)
  return(asv_data)
  
}

infer_asv_command <-function(intermediate_path, denoised_data_path){
  denoised_data_path <- file.path(intermediate_path, denoised_data_path)
  if (file.exists(denoised_data_path)) {
    print("File already exists")
  } else {
    run_dada <- function(direction) {
      lapply(unique(cutadapt_data$primer_name), function(primer_name) infer_asvs(primer_name, direction)) %>%
        unlist(recursive = FALSE)
    }
    dada_forward <- run_dada("Forward")
    dada_reverse <- run_dada("Reverse")
    save(dada_forward, dada_reverse, file = denoised_data_path)
    print("File is now saved")
  }
}

merge_reads_command <- function(intermediate_path, merged_read_data_path){
  # (Hung) Getting user input for min_Overlap and max_Mismatch
  print("Please enter the minimum length of overlap required for merging the forward and reverse reads (Default: 12).")
  min_Overlap <- as.integer(readline())
  print("\n")
  print("Please enter the maximum mismatches allowed in the overlap region (Default: 0).")
  max_Mismatch <- as.integer(readline())
  
  
  merged_read_data_path <- file.path(intermediate_path, merged_read_data_path)
  formatted_ref_dir <- 
  if (file.exists(merged_read_data_path)) {
    merged_reads <- readRDS(merged_read_data_path)
  } else {
    merged_reads <- mergePairs(dadaF = dada_forward,
                               derepF = file.path('intermediate_data/', 'filtered_sequences', names(dada_forward)),
                               dadaR = dada_reverse,
                               derepR = file.path('intermediate_data/', 'filtered_sequences', names(dada_reverse)),
                               minOverlap = min_Overlap,
                               maxMismatch = max_Mismatch,
                               returnRejects = TRUE, 
                               verbose = TRUE)
    saveRDS(merged_reads, file = merged_read_data_path)
    return(merged_reads)
  }
}

#################################################################################################################
#TODO make some sort of plot to come to look at the amount of overlap and percent identity in the overlap region. This
#could be a check to get an idea of how each locus is getting merged

#FROM ZACH need to make some pretty big modifications
##Need to do some renaming to be able to make plot
#non_empty_merged_reads <- merged_reads[map_dbl(merged_reads, nrow) > 0]
#merge_data <- non_empty_merged_reads %>%
#  bind_rows() %>%
#  mutate(file_name = rep(names(non_empty_merged_reads), map_int(non_empty_merged_reads, nrow)),
#         sample_id = gsub(file_name, pattern = '_.+$', replacement = '')) %>%
#  as_tibble()

#metadata <- read_csv(file.path('intermediate_data', 'metadata.csv'))
#merge_data <- left_join(merge_data, metadata, by = 'sample_id')

#merge_data <- select(merge_data, nmatch, nmismatch, nindel, accept)
#merge_data

#merge_data <- mutate(merge_data,
#                     overlap = nmatch + nmismatch,
#                     mismatch = nmismatch + nindel,
#                    identity = (overlap - mismatch) / overlap)

#merge_plot <- merged_data %>%
#  select(locus, mismatch,  accept, overlap) %>%
#  rename('Locus' = locus, 'Mismatches and Indels' = mismatch, 'Merged' = accept, 'Overlap Length' = overlap) %>%
#  gather(key = 'stat', value = 'value', -Locus, -Merged) %>%
#  ggplot(aes(x = value, fill = Merged)) +
#  facet_grid(Locus ~ stat, scales = 'free') +
#  geom_histogram(bins = 30) +
#  scale_fill_viridis_d(begin = 0.8, end = 0.2) +
#  labs(x = '', y = 'ASV count', fill = 'Merged') +
#  theme(panel.grid.major.x = element_blank(),
#        panel.grid.minor = element_blank(),
#       legend.position="bottom") 
#ggsave(merge_plot, filename = 'read_merging.png', path = 'results', width = 8, height = 8)
#merge_plot


#################################################################################################################
#################################################################################################################
## Name: makeSeqtab; make_abund_table
## Purpose: DADA2 functions to make sequence table 
## Inputs: Cutadapt data table and path to filtered sequences
## Returns: ASV abundance table
##################################################################################################################
makeSeqtab<-function(merged_reads){
  raw_seqtab<-makeSequenceTable(merged_reads)
  return(raw_seqtab)
}

#Remove chimeras and short sequencs
make_abund_table<-function(seqtab){
  seqtab.nochim <- removeBimeraDenovo(raw_seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  asv_abund_table <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) >= min_asv_length]
  return(asv_abund_table)
}
#################################################################################################################
## Name: Assign_tax_dada2-TODO-FUNCTIONS NOT MADE-JUST TESTING FUNCTIONALITY USING ZACH'S SCRIPTS
## Purpose: Core functions are from DADA2, and wrapper function adapted from Zach Foster. 
## Inputs: Cutadapt data table and path to filtered sequences
## Returns: Separate ASV abundance matrices for rps10 and ITS
##################################################################################################################
prep_abund_table <-function(cutadapt_data, asv_abund_table, locus){
  rownames(asv_abund_table) <- sub(rownames(asv_abund_table), pattern = ".fastq.gz", replacement = "")
  cutadapt_data$file_id_primer <- paste0(cutadapt_data$file_id, "_", cutadapt_data$primer_name)
  asv_abund_matrix<- asv_abund_table[rownames(asv_abund_table) %in% cutadapt_data$file_id_primer[cutadapt_data$primer_name == locus], ]
  return(asv_abund_matrix)
}


#OPTIONAL FUNCTION IF YOU HAVE MORE THAN ONE LOCUS YOU ARE WORKING WITH. Generalize even more-do a similar step but all possible loci. This may get pretty complicated if you have more than two loci, but at least user could specify which loci? 
#Need to account for the fact that there potentially 2+ loci. ASVs can only be in one locus and not  both. ASVs should either be in one locus or another.If it is in both, assign it to the one with more reads.
#extra steps if you have more than one locus; Very specific to working with ITS and rps10.
separate_abund_table <- function(abund_asv_its, abund_asv_rps10, intermediate_path, separate_abund_path){
  separate_abund_path <- file.path(intermediate_path, separate_abund_path)
  in_both <- colSums(abund_asv_its) != 0 & colSums(abund_asv_rps10) != 0
  assign_to_its <- in_both & colSums(abund_asv_rps10) > colSums(abund_asv_its)
  assign_to_rps10 <- in_both & colSums(abund_asv_rps10) < colSums(abund_asv_its)
  is_its <- (colSums(abund_asv_its) != 0 & colSums(abund_asv_rps10) == 0) | assign_to_its
  is_rps10 <- (colSums(abund_asv_rps10) != 0 & colSums(abund_asv_its) == 0) | assign_to_rps10
  abund_asv_its <- abund_asv_its[ , is_its]
  abund_asv_rps10 <- abund_asv_rps10[ , is_rps10]
  #The number of ASVs left in the two groups should sum to the total number of ASVs, since there should be no overlap.
  stopifnot(ncol(abund_asv_its) + ncol(abund_asv_rps10) == ncol(asv_abund_table)) #make more messaging
  save(abund_asv_its, abund_asv_rps10, file = separate_abund_path)
}

#################################################################################################################
## Name: CORE DADA2 FUNCTION-Assign taxonomy
## Purpose: Core functions are from DADA2
## Inputs: asv abundance table from previous step
## Returns: Taxonomic table for each locus
##################################################################################################################

##THESE DADA2 FUNCTIONS CAN BE SLOW-RPS10 database is much smaller, but ITS is quite large. 
##On my personal computer, I actually had to use a smaller version of the ITS database.
##generalize more

assign_taxonomyDada2<-function(abund_asv_table, ref_database){
  tax_results<- assignTaxonomy(abund_asv_table, 
                               refFasta = file.path(intermediate_path, 'reference_databases', ref_database), 
                               taxLevels = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Reference"),
                               minBoot = 0,
                               tryRC = TRUE,
                               outputBootstraps = TRUE,
                               multithread = TRUE)
  return(tax_results)
}

#################################################################################################################
## Align to reference sequence for percent ID. Functions are entirely from Zach
## Purpose: As a quality check, align ASV sequences to reference sequence to assign a percent identity. A high bootstrap value doesn't always  mean a good match, just that it is better than other matches. 
## Inputs: Input database, taxonomic results from above
## Returns: Percent identity of ASV sequences relative to reference sequence. Final function will add PIDs into taxonomic assignment results as another rank. You will then see percent identity of ASV to reference sequence. 
##################################################################################################################

get_ref_seq <- function(tax_result, db) {
  ref_i <- as.integer(str_match(tax_result$tax[, 'Reference'], '^.+_([0-9]+)$')[ ,2])
  db[ref_i]
}

get_align_pid <- function(ref, asv) {
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  align <-  pairwiseAlignment(pattern = asv, subject = ref, type = 'global-local')
  is_match <- strsplit(as.character(align@pattern), '')[[1]] == strsplit(as.character(align@subject), '')[[1]]
  sum(is_match) / length(is_match)
}

get_pids <- function(tax_result, db) {
  ref_seq <- get_ref_seq(tax_result, db)
  asv_seq <- rownames(tax_result$tax)
  future_map2_dbl(ref_seq, asv_seq, get_align_pid) * 100
}

add_pid_to_tax <- function(tax_result, pid) {
  tax_result$tax <- cbind(tax_result$tax, ASV = rownames(tax_result$tax))   
  tax_result$boot <- cbind(tax_result$boot, ASV = pid)
  return(tax_result)
}

#I will combine the taxonomic assignments and bootstrap values for each locus into a single classification vector.
#This will store all the taxonomic and bootstrap information in a single vector. 
#################################################################################################################
## Combine taxonomic assignments and bootstrap values for each locus into single vector. 
## Purpose: For downstream work, have taxonomic classifications for each locus in a single vector
## Inputs: Each of the taxonomic assignments and bootstrap values for each locus. 
## Returns: All taxonomic and bootstrap info is in a single vector
##################################################################################################################
assignTax_as_char <- function(tax_result) {
  out <- vapply(1:nrow(tax_result$tax), FUN.VALUE = character(1), function(i) {
    paste(tax_result$tax[i, ],
          tax_result$boot[i, ],
          colnames(tax_result$tax), 
          sep = '--', collapse = ';')
  })
  names(out) <- rownames(tax_result$tax)
  return(out)
  #check
  stopifnot(all(names(out) %in% colnames(asv_abund_table)))
  stopifnot(all(! duplicated(names(out))))
}


#Reformat ASV table
format_abund_matrix <- function(asv_abund_table) {
  formatted_abund_asv <- t(asv_abund_table)
  colnames(formatted_abund_asv) <- sub(colnames(formatted_abund_asv), pattern = ".fastq.gz$", replacement = "")
  formatted_abund_asv <- cbind(sequence = rownames(formatted_abund_asv), 
                              taxonomy = seq_tax_asv[rownames(formatted_abund_asv)], 
                              formatted_abund_asv)
  formatted_abund_asv <- as_tibble(formatted_abund_asv)
  write_csv(formatted_abund_asv, file = file.path(intermediate_path, 'final_asv_abundance_table.csv'))
  #make results folder
  print(formatted_abund_asv)
  return(formatted_abund_asv)
}

#################################################################################################################
## Name: TODO Take inventory of read counts after various steps and make plots to visualize the outputs along the way
## Purpose: TODO
## Inputs: TODO
## Returns: TODO

#Zach has one example of this, but I think we an simplify this. 
##################################################################################################################

#################################################################################################################
## Name: Convert to Phyloseq object and Taxa objects for downstream work
## Purpose: TODO
## Inputs: TODO
## Returns: TODO
#################################################################################################################

#****************************************************************FUNCTION CALLS************************************************************************
#In progress
#Still trying to figure out best ways to group functions since there are MANY small functions 
#this is all for testing - TODO eventually create/edit a main function
#primer_data-just need primer_name, forward, reverse
#metadata-just need sample (needs to be first), primer name
raw_path <-"fastq_reads/Trimmed_function_reads/"
primer_path <-file.path(raw_path, "raw_data/primer_info_example.csv")
metadata_path <-file.path(raw_path, "raw_data/metadata_example.csv")
#fastq_path<-"~/rhod_metabarcoding/raw_data" just raw data path
# to large to fit on GitHub - example fastq files on local machine

#create intermediate in working directory
primer_data <- prepare_primers(primer_path)
metadata <- prepare_metadata(metadata_path, primer_data)

create_intermediate("intermediate_data/") 
intermediate_path <- "intermediate_data/"


#this may take while

fastq_data <- prepare_fastq(raw_path,intermediate_path, metadata)

#this also could take a hot sec
pre_primer_hit_data<- pre_primer_hit_data(primer_data, fastq_data, intermediate_path)
#be careful about unzipped files. Output was bad
print(pre_primer_hit_data)

pre_primer_plot <- primer_hit_plot(pre_primer_hit_data, fastq_data, metadata)
print(pre_primer_plot)
ggsave(pre_primer_plot, filename = 'pre_trim_counts.pdf', path = intermediate_path, width = 8, height = 8)

cutadapt_data <- cutadapt_tibble(fastq_data, metadata, intermediate_path)

#running the actual cutadapt program
cutadapt_run("/Users/phanhung2/.local/bin/cutadapt", cutadapt_data)

#incorporate other quality plots? 
#filter steps

#check qc before trim
quality_plots<-plot_qc(cutadapt_data, "S1")
ggsave(quality_plots, filename = 'Prefilter_qc.pdf', path = intermediate_path, width = 8, height = 8) #incorporate into fxn

filter_and_trim(intermediate_path, cutadapt_data)

#checking if cutadapt works (this may take awhile)
post_primer_hit_data <- post_primer_hit_data(primer_data, cutadapt_data, intermediate_path) #incorporate into fxns
print(post_primer_hit_data)
post_primer_plot <- primer_hit_plot(post_primer_hit_data, new_fastq_data, metadata)
print(post_primer_plot)
ggsave(post_primer_plot, filename = 'Post_trim_counts.pdf', path = intermediate_path, width = 8, height = 8) #incorporate into fxn

quality_plots2<-post_trim_qc(cutadapt_data, "S1")
ggsave(quality_plots2, filename = 'Filtered_qc.pdf', path = intermediate_path, width = 8, height = 8) #incorporate into fxn

#database functions
##Make sure you you have fasta files of your databases first
#first manually move your downloaded database files to the reference databases folder in intermediate data directory-could automate this or retrieve from downloads folder? 

create_ref_database(intermediate_path)
rps10_data<-format_database(raw_path, "oomycetedb.fasta")
its_data<-format_database2(raw_path, "unite.fasta")

#####NOT FINISHED YET######
#ASV commands-still not that generalized-initial proof of concept-not complete
#This part is a bit clunky to run command and then load saved files. Is there a better solution. Something like Snakemake or Nextflow solves some problems, but not compatible with R? 
infer_asv_command(intermediate_path, "Denoised_data.Rdata")
denoised_data_path <- file.path(intermediate_path, "Denoised_data.Rdata")
load(denoised_data_path)


#merge reads
merged_reads<-merge_reads_command("~/rhod_metabarcoding4/intermediate_data", "Merged_reads.rds")
#Get more info on merged reads

#Sequence tables
raw_seqtab<-makeSeqtab(merged_reads) #incorporate into functions
dim(raw_seqtab)
hist(nchar(getSequences(raw_seqtab)))

asv_abund_table<-make_abund_table(raw_seqtab)
print(sum(asv_abund_table)/sum(raw_seqtab)) #incorporate into functions
dim(asv_abund_table)
hist(nchar(getSequences(asv_abund_table)))

#make abundance table for each locus present within each sequence
rps10_abund_asv <- prep_abund_table(cutadapt_data, asv_abund_table, "rps10")
its_abund_asv <- prep_abund_table(cutadapt_data, asv_abund_table, "ITS")

#Ensure that abundance tables are separate for each locus
separate_abund_table(its_abund_asv, rps10_abund_asv, "~/rhod_metabarcoding4/intermediate_data", "Separate_abund.Rdata")
separate_abund_path <- file.path("~/rhod_metabarcoding4/intermediate_data/", "Separate_abund.Rdata")
load(separate_abund_path)

#Assign taxonomy for each locus using DADA2 function
#need to fix
tax_results_rps10_asv <- assign_taxonomyDada2(abund_asv_rps10, "rps10_reference_db.fa")
#Note for ITS ref database, I refer to a much smaller version of the reformatted ITS database because my computer couldn't handle the larger one. I've included this in the Google drive folder. After the new reference_databases directory is created above, I just placed the its_short.fasta into this directory too. You computer may not have an issue.
tax_results_its_asv <- assign_taxonomyDada2(abund_asv_its, "its_short.fasta")


#Percent identity as quality check and for potential downstream filtering
#need to generalize
its_seqs <- read_fasta(file.path(intermediate_path, "reference_databases",  'its_short.fasta'))
rps10_seqs <- read_fasta(file.path(intermediate_path, "reference_databases", 'rps10_reference_db.fa'))
rps10_pids_asv <- get_pids(tax_results_rps10_asv, rps10_seqs)
its_pids_asv <- get_pids(tax_results_its_asv, its_seqs)


tax_results_rps10_asv <- add_pid_to_tax(tax_results_rps10_asv, rps10_pids_asv)
tax_results_its_asv <- add_pid_to_tax(tax_results_its_asv, its_pids_asv)

seq_tax_asv <- c(assignTax_as_char(tax_results_rps10_asv), assignTax_as_char(tax_results_its_asv))
head(seq_tax_asv)

formatted_abund_asv<-format_abund_matrix(asv_abund_table)

#final function to combine both results
#plot to visualize pre and post filtering results
#taxa and phyloseq objects for downstream work? (would require certain minimum amount of metadata)

#****************************************************************NOTES************************************************************************
##To consider-look into any other informative plots as quality control checks
#Incude figure that incorporates read count distributions by length
#Discuss much control should user have over parameters?
#Review Cutadapt parameters
#BLAST taxonomic comparison?
#Save output plots in a results directory and have a message where this output file is in case people don't use Rstudio***
#Make sure paths are consistent 
#what functions still need to be added? 
#What best practices are violated? 