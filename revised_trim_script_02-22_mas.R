#load packages 
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
library(sessioninfo)

seed <- 1
set.seed(seed)
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

#################################################################################################################
## Name: prepare_metadata
## Purpose: to read in the metadata from user and combine it with the primer data
## Inputs: the metadata path, primer data tibble created in prepare_primers()
## Returns: a data frame containing the metadata and primer data now combined
##################################################################################################################
prepare_metadata <- function(metadata_path, primer_data){
  metadata <- read_csv(metadata_path) %>%
    left_join(primer_data, by = c("primer_name"))
  #this could cause issues in other data sets, if the name of the sample is no the first column
  names(metadata)[1] <- "sample_id"
  return(metadata)
}

#################################################################################################################
## Name: create_intermediate
## Purpose: to create a file in the current working directory used for program files later on
## Inputs: none
## Returns: none 
##################################################################################################################
create_intermediate <- function(){
  if (! dir.exists("intermediate_data")) {
    dir.create("intermediate_data")
  }
}

#Testing function
#intermediate <-create_intermediate("~/rhod_metabarcoding4")
#################################################################################################################
## Name: read_fastq 
## Purpose: take in the fastq files from the user and creating a tibble that organizes them
## Inputs: fastq file path 
## Returns: a tibble with the fastq file paths, the direction of the sequences, and names of sequences
##################################################################################################################
read_fastq <- function(fastq_path){
  fastq_paths <- list.files(fastq_path, pattern = "\\.fastq")
  #constructing a tibble - special data frame with improved behaviors. 
  fastq_data <- tibble(file_id = sub(fastq_paths, pattern = "\\.fastq\\.gz$", replacement = ""),
                       sample_id = sub(fastq_paths, pattern = "_.+$", replacement = ""),
                       #direction: the grep1 is looking to match the pattern, relates a TRUE or FALSE with it (aka 1 or 0) and adds 1 to find the 
                       #correct index to put Reverse or Forward
                       direction = c( "Reverse", "Forward")[BiocGenerics::grepl(fastq_paths, pattern = "_R1") + 1],
                       raw_path = file.path(fastq_path, fastq_paths))
  
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
#################################################################################################################
## Name: prepare_fastq (don't run above function individually)
## Purpose: does all of the processing/editing of the fastqs in order to make ready for cutadapt
## Inputs: path to the fastq files, the intermediate data path
## Outputs: a data frame with the fastq information 
##################################################################################################################
prepare_fastq <- function(fastq_path,intermediate_path, metadata){
  fastq_data <- read_fastq(fastq_path)
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
}

#################################################################################################################
## Name: cutadapt_tibble
## Purpose: create a tibble that can be fed into cutadapt 
## Inputs: a fastq_data tibble created by prepare_fastq
## Returns: a cutadapt tibble
##################################################################################################################
cutadapt_tibble <- function(fastq_data, metadata, intermediate_path){ #new_fastq_data needed, why not just fastq_data
  cutadapt_data <- metadata %>%
    left_join(fastq_data, by = c("sample_id" = "sample_id"))
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
#Not exactly like Sam's function. Needs additional simplification and generalization
cutadapt_run <- function(cutadapt_path, cutadapt_data){
    cutadapt <- cutadapt_path
    tryCatch(system2(cutadapt, args = "--version"),
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
  
  #simplify like Sam
  #we are jsut building a column in cutadapt data that has the command line arguements in it
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
  command_args= paste(R1_flags, R2_flags, "-n", 2, "-o",  fwd_trim, "-p", rev_trim, "--minimum-length", 50,   "--untrimmed-output", fwd_untrim, "--untrimmed-paired-output", rev_untrim,"--quiet",fwd_prefilt, rev_prefilt)
  if (! all(file.exists(c(cutadapt_data$trimmed_path,cutadapt_data$untrimmed_path)))) {
    cutadapt_output <- future_map(command_args, ~system2(cutadapt, args = .x))
  }
  
  return(cutadapt_data)
  
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

#################################################################################################################
## Name: filter_and_trim
## Purpose: to filter out reads with Ns, below a length of 50, or greater than expected error limis. References filterAndTrim function from DADA2.
## Inputs: intermediate file path, cutadapt tibble
## Returns: none
##################################################################################################################
filter_and_trim <- function(intermediate_path, cutadapt_data){
  #make a place to put the results 
  expected_error_filter_limit <- 5 #default setting
  truncation_qual_limit <- 5 #default setting 
  filtered_read_dir <- file.path(intermediate_path, "filtered_sequences")
  cutadapt_data$filtered_path <- file.path(filtered_read_dir, paste0(cutadapt_data$file_id, "_", cutadapt_data$primer_name, ".fastq.gz"))
  if(! all(file.exists(cutadapt_data$filtered_path))){
    filter_results <- filterAndTrim(fwd = cutadapt_data$trimmed_path[cutadapt_data$direction == "Forward"], 
                                    filt = cutadapt_data$filtered_path[cutadapt_data$direction == "Forward"],
                                    rev =  cutadapt_data$trimmed_path[cutadapt_data$direction == "Reverse"], 
                                    filt.rev = cutadapt_data$filtered_path[cutadapt_data$direction == "Reverse"], 
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

#better visualized as a plot-merge different ? 
#################################################################################################################
## Name: post_trim_qc
## Purpose: Function is from DADA2
## Inputs: Sample you are interested in plotting (use quotes), cutadapt tibble
## Returns: Plots of quality score at each base position after trimming
##################################################################################################################
post_trim_qc<-function(cutadapt_data, pick_sample){
  #just retrieve all plots for first sample
  sample_info = cutadapt_data$filtered_path[cutadapt_data$sample_id == pick_sample]
  quality_plots<-plotQualityProfile(sample_info)
  return(quality_plots)
}

#****************************************************************FUNCTION CALLS************************************************************************
#In progress
#From Martha-still trying to figure out best ways to run sequences of functions
#this is all for testing - TODO eventually create/edit a main function
primer_path <- "~/rhod_metabarcoding4/primer_info_m.csv"
metadata_path <- "~/rhod_metabarcoding4/metadata.csv"
# to large to fit on GitHub - example fastq files on local machine
fastq_path <- "~/rhod_metabarcoding4"
create_intermediate()
intermediate_path <- ("~/rhod_metabarcoding4/intermediate_data")
create_intermediate()

primer_data <- prepare_primers(primer_path)
metadata <- prepare_metadata(metadata_path, primer_data)
#this may take while
fastq_data <- prepare_fastq(fastq_path, intermediate_path)
#this also could take a hot sec
pre_primer_hit_data <- pre_primer_hit_data(primer_data, fastq_data, intermediate_path)
print(pre_primer_hit_data)

pre_primer_plot <- primer_hit_plot(pre_primer_hit_data, fastq_data, metadata)
print(pre_primer_plot)

cutadapt_data <- cutadapt_tibble(fastq_data, metadata, intermediate_path)

#running the actual cutadapt program
cutadapt_run("/Users/masudermann/.local/bin/cutadapt", cutadapt_data)

#filter steps
#cutadapt_data1 <-filter_and_trim("~/rhod_metabarcoding4/intermediate_data", cutadapt_data1)
cutadapt_data <- filter_and_trim(intermediate_path, cutadapt_data)

#checking if cutadapt works (this may take awhile)
post_primer_hit_data <- post_primer_hit_data(primer_data, cutadapt_data, intermediate_path)
print(post_primer_hit_data)
post_primer_plot <- primer_hit_plot(post_primer_hit_data, new_fastq_data, metadata)
print(post_primer_plot)

post_trim_qc(cutadapt_data, "S1")
##To consider
#any other informative plots
#After discussion with Ricardo-Read count distributions by length
#Better assess read quality before filter steps
