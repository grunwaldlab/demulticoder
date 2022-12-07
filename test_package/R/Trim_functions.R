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
make_cutadapt_tibble <- function(fastq_data, metadata, intermediate_path){ #new_fastq_data needed, why not just fastq_data
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

#rename min_length
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
run_cutadapt <- function(cutadapt_path, cutadapt_data, min_length=20){
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
  command_args= paste(R1_flags, R2_flags, "-n", 2,  "-o",  fwd_trim, "-p", rev_trim, "--minimum-length", min_length,"--untrimmed-output", fwd_untrim, "--untrimmed-paired-output", rev_untrim,"--quiet", fwd_prefilt, rev_prefilt)
  if (! all(file.exists(c(cutadapt_data$trimmed_path)))) {
    cutadapt_output <- furrr::future_map(command_args, ~system2(cutadapt, args = .x))
  }

}
