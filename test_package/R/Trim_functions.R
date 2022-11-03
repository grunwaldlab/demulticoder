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
primer_hit_plot <- function(primer_hits, fastq_data, intermediate_path, plot_name){
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

#' Make a barplot of primers identified on reads after trim steps
#'
#' @param primer_hits A number of reads in which the primer is found
#' @param fastq_data A tibble with the fastq file paths, the direction of the sequences, and names of sequences
#' @param metadata A metadata containing the concatenated metadata and primer data
#'
#' @return
#' @export
#'
#' @examples
post_primer_hit_plot <- function(primer_hits, fastq_data, intermediate_path, plot_name){
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

#' Prepare for primmer trimming with Cutaapt. Make new sub-directories and specify paths for the trimmed and untrimmed reads
#'
#' @param fastq_data A path to FASTQ files for analysis
#' @param metadata A metadata containing the concatenated metadata and primer data
#' @param intermediate_path A path to the intermediate folder and directory
#'
#' @return
#' @export
#'
#' @examples
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

#' Core function for running cutadapt
#'
#' @param cutadapt_path A path to the cutadapt program
#' @param cutadapt_data Intermediate_data folder with trimmed and filtered reads for each sample
#' @param min_length Read lengths that are lower than this threshold will be discarded. Default is 50.
#'
#' @return
#' @export
#'
#' @examples
cutadapt_run <- function(cutadapt_path, cutadapt_data, min_length=50){
  cutadapt <- cutadapt_path
  tryCatch(system2(cutadapt, args = "--version"),
           warning=function(w){
             stop("cutadapt cannot be found on PATH. Check if it is installed?")
           })
  #simplify
  #may need to demultiplex first-meed to assess this
  cutadapt <- path.expand(cutadapt_path)
  min_length=min_length
  R1_flags = unique(paste("-g",cutadapt_data$forward, "-a", cutadapt_data$r_rc))
  R2_flags = unique(paste("-G",cutadapt_data$reverse, "-A", cutadapt_data$f_rc))
  fwd_trim = cutadapt_data[cutadapt_data$direction == "Forward", ][["trimmed_path"]]
  rev_trim = cutadapt_data[cutadapt_data$direction == "Reverse", ][["trimmed_path"]]
  fwd_untrim = cutadapt_data[cutadapt_data$direction == "Forward", ][["untrimmed_path"]]
  rev_untrim = cutadapt_data[cutadapt_data$direction == "Reverse", ][["untrimmed_path"]]
  fwd_prefilt = cutadapt_data[cutadapt_data$direction == "Forward", ][["prefiltered_path"]]
  rev_prefilt = cutadapt_data[cutadapt_data$direction == "Reverse", ][["prefiltered_path"]]
  #consider parameters
  command_args= paste(R1_flags, R2_flags, "-n", 2,  "-o",  fwd_trim, "-p", rev_trim, "--minimum-length", min_length,"--minimum-length", min_length,"--untrimmed-output", fwd_untrim, "--untrimmed-paired-output", rev_untrim,"--quiet", fwd_prefilt, rev_prefilt)
  if (! all(file.exists(c(cutadapt_data$trimmed_path)))) {
    cutadapt_output <- furrr::future_map(command_args, ~system2(cutadapt, args = .x))
  }

}
