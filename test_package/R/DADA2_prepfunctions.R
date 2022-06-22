#' Wrapper function for core DADA2 filter and trim function for first filtering step
#'
#' @param fastq_data
#' @param intermediate_path
#' @param metadata
#' @inheritParams filterAndTrim
#'
#' @return
#' @export
#'
#' @examples
remove_ns <- function(fastq_data, intermediate_path){
  prefiltered_read_dir <- file.path(intermediate_path, "prefiltered_sequences")
  fastq_data$prefiltered_path <- file.path(prefiltered_read_dir, base::basename(fastq_data$raw_data_path))
  #if the files do not exist in the prefiltered path path (which clearly they don't)
  #raw path is the path to the actual fastq files in mock community
  #fwd takes from the raw path data
  #filt is the path to the output filtered files that we created
  #the sequences it chose to put into prefiltered path had to have no Ns in them
  if(! all(file.exists(fastq_data$prefiltered_path))){
    dada2::filterAndTrim(fwd = fastq_data[fastq_data$direction == "Forward", ][["raw_data_path"]],
                         filt = fastq_data[fastq_data$direction == "Forward", ][["prefiltered_path"]],
                         rev = fastq_data[fastq_data$direction == "Reverse", ][["raw_data_path"]],
                         filt.rev = fastq_data[fastq_data$direction == "Reverse", ][["prefiltered_path"]],
                         maxN = 0,
                         multithread = TRUE)
  }
  return(fastq_data)


}

#' Wrapper function for plotQualityProfile function
#'
#' @param cutadapt_data
#' @param intermediate_path
#' @inheritParams plotQualityProfile
#'
#' @return
#' @export
#'
#' @examples
plot_qc<-function(cutadapt_data, intermediate_path){
  #just retrieve all plots for first sample
  for (i in unique(cutadapt_data$sample_id))
  {
    sample_info = cutadapt_data$trimmed_path[cutadapt_data$sample_id == i]
    quality_plots<-dada2::plotQualityProfile(sample_info)
    name1=paste0('qcpre_trim_plot_', i,'.pdf')
    ggplot2::ggsave(quality_plots, filename = name1, path = intermediate_path, width = 8, height = 8)
    #print(quality_plots)
  }
}


#' Wrapper function for filterAndTrim function from DADA2 after primer removal
#'
#' @param intermediate_path
#' @param cutadapt_data
#'
#' @return
#' @export
#'
#' @examples
filter_and_trim <- function(intermediate_path, cutadapt_data){
  expected_error_filter_limit <- 5 #default setting
  truncation_qual_limit <- 5 #default setting
  min_Length <- 50
  max_Length <- 285
  filtered_read_dir <- file.path(intermediate_path, "filtered_sequences")
  cutadapt_data$filtered_path <- file.path(filtered_read_dir, paste0(cutadapt_data$file_id, "_", cutadapt_data$primer_name, ".fastq.gz"))
  if(! all(file.exists(cutadapt_data$filtered_path))){
    filter_results <- dada2::filterAndTrim(fwd = cutadapt_data$trimmed_path[cutadapt_data$direction == "Forward"],
                                           filt = cutadapt_data$filtered_path[cutadapt_data$direction == "Forward"],
                                           rev =  cutadapt_data$trimmed_path[cutadapt_data$direction == "Reverse"],
                                           filt.rev = cutadapt_data$filtered_path[cutadapt_data$direction == "Reverse"],
                                           maxN = 0,
                                           maxEE = c(expected_error_filter_limit, expected_error_filter_limit),
                                           truncQ = truncation_qual_limit,
                                           minLen = min_Length,
                                           maxLen = max_Length,
                                           rm.phix = TRUE,
                                           compress = TRUE,
                                           matchIDs = TRUE,
                                           multithread = TRUE)
    filter_results <- as_tibble(filter_results)
    filter_results_path<-file.path(intermediate_path, "filter_results.RDS")
    saveRDS(filter_results, file = filter_results_path)
    filter_results_path2<-file.path(intermediate_path, "filter_results.RData")
    save(filter_results, file = filter_results_path2)
    print(colMeans(filter_results))
  }
}


#' Wrapper script for plotQualityProfile after trim steps and primer removal.
#'
#' @param cutadapt_data
#' @param intermediate_path
#'
#' @return
#' @export
#'
#' @examples
post_trim_qc<-function(cutadapt_data, intermediate_path){
  #just retrieve all plots for first sample
  for (i in unique(cutadapt_data$sample_id))
  {
    sample_info = cutadapt_data$trimmed_path[cutadapt_data$sample_id == i]
    quality_plots2<-dada2::plotQualityProfile(sample_info)
    name=paste0('qcpost_trim_plot_', i,'.pdf')
    ggplot2::ggsave(quality_plots2, filename = name, path = intermediate_path, width = 8, height = 8)
    #print(quality_plots2)
  }
}
