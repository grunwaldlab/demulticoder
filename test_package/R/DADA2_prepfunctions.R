#' Wrapper function for core DADA2 filter and trim function for first filtering step
#' @inheritParams dada2::filterAndTrim
#' @param fastq_data A tibble with the fastq file paths, the direction of the sequences, and names of sequences
#' @param intermediate_path A path to the intermediate folder
#' @param metadata A metadata containing the concatenated metadata and primer data
#' @inheritParams filterAndTrim
#'
#' @return
#' @export
#'
#' @examples
remove_ns <- function(fastq_data, intermediate_path, maxN= 0, multithread = TRUE){
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
                         maxN=maxN,
                         multithread=multithread)
  }
  return(fastq_data)
}

#' Wrapper function for plotQualityProfile function
#' @inheritParams dada2::plotQualityProfile
#' @param cutadapt_data Intermediate_data folder with trimmed and filtered reads for each sample
#' @param intermediate_path A path to the intermediate folder
#' @inheritParams plotQualityProfile
#'
#' @return
#' @export
#'
#' @examples
plot_qc<-function(cutadapt_data, intermediate_path, n=500000){
  #just retrieve all plots for first sample
  for (i in unique(cutadapt_data$sample_id))
  {
    sample_info = cutadapt_data$trimmed_path[cutadapt_data$sample_id == i]
    quality_plots<-dada2::plotQualityProfile(sample_info, n)
    name1=paste0('qcpre_trim_plot_', i,'.pdf')
    ggplot2::ggsave(quality_plots, filename = name1, path = intermediate_path, width = 8, height = 8)
    #print(quality_plots)
  }
}


#' Wrapper function for filterAndTrim function from DADA2 after primer removal
#' @inheritParams dada2::filterAndTrim
#' @param intermediate_path A path to the intermediate folder
#' @param cutadapt_data Intermediate_data folder with trimmed and filtered reads for each sample
#'
#' @return
#' @export
#'
#' @examples
filter_and_trim <- function(intermediate_path, cutadapt_data,  maxEE = Inf, truncQ = 2, minLength = 20, maxLength = Inf, truncLen = 0, maxN = 0, minQ=0, rm.phix=TRUE, multithread=FALSE, matchIDs=FALSE, verbose=FALSE, qualityType="Auto", OMP=TRUE, n=1e+05,id.sep="\\s", rm.lowcomplex=0, orient.fwd=NULL, id.field=NULL){
  filtered_read_dir <- file.path(intermediate_path, "filtered_sequences")
  cutadapt_data$filtered_path <- file.path(filtered_read_dir, paste0(cutadapt_data$file_id, "_", cutadapt_data$primer_name, ".fastq.gz"))
  if(! all(file.exists(cutadapt_data$filtered_path))){
    filter_results <- dada2::filterAndTrim(fwd = cutadapt_data$trimmed_path[cutadapt_data$direction == "Forward"],
                                         filt = cutadapt_data$filtered_path[cutadapt_data$direction == "Forward"],
                                         rev =  cutadapt_data$trimmed_path[cutadapt_data$direction == "Reverse"],
                                         filt.rev = cutadapt_data$filtered_path[cutadapt_data$direction == "Reverse"],
                                         maxN = maxN,
                                         maxEE = c(maxEE, maxEE),
                                         truncLen=truncLen,
                                         truncQ = truncQ,
                                         minLen = minLength,
                                         maxLen = maxLength,
                                         minQ = minQ,
                                         trimLeft=0,
                                         trimRight=0,
                                         rm.phix = rm.phix,
                                         compress = TRUE,
                                         matchIDs = FALSE,
                                         multithread = multithread,
                                         verbose=verbose,
                                         OMP=TRUE,
                                         n=1e+05,
                                         id.sep="\\s",
                                         rm.lowcomplex=0,
                                         orient.fwd=NULL,
                                         qualityType="Auto",
                                         id.field=NULL)
    filter_results <- as_tibble(filter_results)
    filter_results_path<-file.path(intermediate_path, "filter_results.RDS")
    saveRDS(filter_results, file = filter_results_path)
    filter_results_path2<-file.path(intermediate_path, "filter_results.RData")
    save(filter_results, file = filter_results_path2)
    print(colMeans(filter_results))
  }
}


#' Wrapper script for plotQualityProfile after trim steps and primer removal.
#' @inheritParams dada2::plotQualityProfile
#' @param cutadapt_data Intermediate_data folder with trimmed and filtered reads for each sample
#' @param intermediate_path A path to the intermediate folder
#'
#' @return
#' @export
#'
#' @examples
post_trim_qc<-function(cutadapt_data, intermediate_path, n=500000){
  #just retrieve all plots for first sample
  for (i in unique(cutadapt_data$sample_id))
  {
    sample_info2 = cutadapt_data$filtered_path[cutadapt_data$sample_id == i]
    quality_plots2<-dada2::plotQualityProfile(sample_info2, n)
    name=paste0('qcpost_trim_plot_', i,'.pdf')
    ggplot2::ggsave(quality_plots2, filename = name, path = intermediate_path, width = 8, height = 8)
    #print(quality_plots2)
  }
}
