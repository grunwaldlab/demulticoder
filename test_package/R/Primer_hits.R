#' Get primer counts fo reach sample before primer removal and trimming steps
#' TODO: don't make if already
#'
#' @param primer_data
#' @param fastq_data
#' @param intermediate_path
#'
#' @return
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
#' @param primer_data
#' @param cutadapt_data
#' @param intermediate_path
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
#' @param primer_hits
#' @param fastq_data
#' @param plot_name
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
#' @param primer_hits
#' @param fastq_data
#' @param metadata
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
