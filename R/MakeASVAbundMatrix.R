#' Make an amplified sequence variants abundance matrix with data processed through preceding steps
#'
#' @param rawSeqTab_fileName A filename as which the raw sequence table will be saved
#' @param abundMatrix_fileName A filename as which the abundance matrix will be saved
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primer_info files
#' @param directory_path_temp User-defined temporary directory to place reads throughout the workflow
#' metadata, and primer_info files
#' @inheritParams infer_asv_command
#' @inheritParams merge_reads_command
#' @return The asv abundance matrix asv_abund_matrix
#' @export make_asv_abund_matrix
#' @examples
#' directory_path<-"~/rps10package/raw_data/rps10_ITS"
#' primer_path <-file.path(directory_path, "primer_info.csv")
#' metadata_path <-file.path(directory_path,"metadata.csv")
#' cutadapt_path<-"/opt/homebrew/bin/cutadapt"
#' data_tables <-
#' prepare_reads(
#' directory_path,
#' primer_path,
#' metadata_path,
#' maxN = 0,
#' )
#' cut_trim(
#' directory_path,
#' cutadapt_path,
#' verbose = TRUE,
#' maxEE = 2,
#' truncQ = 5,
#' minLen = 200,
#' maxLen = 297,
#' minCutadaptlength = 50
#')
#' asv_abund_matrix <-
#' make_asv_abund_matrix(
#' directory_path,
#' minOverlap = 15,
#' maxMismatch = 2,
#' verbose = TRUE,
#' )
#'
#'
#'
#'
make_asv_abund_matrix <- function(directory_path,
                                  multithread = FALSE,
                                  nbases = 1e+08,
                                  errorEstimationFunction = loessErrfun,
                                  randomize = FALSE,
                                  MAX_CONSIST = 10,
                                  OMEGA_C = 0,
                                  qualityType = "Auto",
                                  nominalQ = FALSE,
                                  obs = TRUE,
                                  err_out = TRUE,
                                  err_in = FALSE,
                                  pool = FALSE,
                                  selfConsist = FALSE,
                                  verbose = FALSE,
                                  minOverlap = 12,
                                  maxMismatch = 0,
                                  returnRejects = FALSE,
                                  justConcatenate = FALSE,
                                  trimOverhang = FALSE,
                                  orderBy = "abundance",
                                  method = "consensus",
                                  min_asv_length = 50) {
  infer_asv_command(
    directory_path,
    directory_path_temp,
    multithread = multithread,
    nbases = nbases,
    errorEstimationFunction = errorEstimationFunction,
    randomize = randomize,
    MAX_CONSIST = MAX_CONSIST,
    OMEGA_C = OMEGA_C,
    qualityType = qualityType,
    nominalQ = nominalQ,
    obs = obs,
    err_out = err_out,
    err_in = err_in,
    pool = pool,
    selfConsist = selfConsist,
    verbose = verbose
  )
  merged_reads <-
    merge_reads_command(
      directory_path,
      minOverlap = minOverlap,
      maxMismatch = maxMismatch,
      returnRejects = returnRejects,
      justConcatenate = justConcatenate,
      verbose = verbose
    )
  countOverlap(merged_reads, directory_path)
  raw_seqtab <- makeSeqtab(merged_reads, orderBy = orderBy)
  asv_abund_matrix <-
    make_abund_matrix(raw_seqtab, min_asv_length = min_asv_length)
  make_seqhist(asv_abund_matrix)
  return(asv_abund_matrix)
}

#' Retrieve the paths of the filtered and trimmed Fastq ifiles
#'
#' @param my_direction Whether primer is in forward or reverse direction
#' @param my_primer_pair_id The the specific barcode id 
#' @param cutadapt_data directory_data folder with trimmed and filtered reads for each sample
#' @keywords internal
get_fastq_paths <- function(my_direction, my_primer_pair_id) {
  data_tables$cutadapt_data %>%
    filter(
      direction == my_direction,
      primer_name == my_primer_pair_id,
      file.exists(filtered_path)
    ) %>%
    pull(filtered_path)
}

#' Core DADA2 function to learn errors and infer ASVs
#'
#' @inheritParams dada2::learnErrors
#' @inheritParams dada2::dada
#' @inheritParams dada2::plotErrors
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primer_info files
#' @param my_primer_pair_id The the specific barcode id 
#' @param my_direction Location of read files and metadata file
#' @return asv_data
#' @keywords internal
infer_asvs <-
  function(directory_path,
           my_primer_pair_id,
           my_direction,
           multithread = FALSE,
           nbases = 1e+08,
           errorEstimationFunction = loessErrfun,
           randomize = FALSE,
           MAX_CONSIST = 10,
           OMEGA_C = 0,
           qualityType = "Auto",
           nominalQ = FALSE,
           obs = TRUE,
           err_out = TRUE,
           err_in = FALSE,
           pool = FALSE,
           selfConsist = FALSE,
           verbose = FALSE) {
    #may need to adjust the parameters
    fastq_paths <- get_fastq_paths(my_direction, my_primer_pair_id)
    error_profile <-
      dada2::learnErrors(
        fastq_paths,
        multithread = multithread,
        nbases = nbases,
        errorEstimationFunction = errorEstimationFunction,
        randomize = randomize,
        MAX_CONSIST = MAX_CONSIST,
        OMEGA_C = OMEGA_C,
        qualityType = qualityType,
        verbose = verbose
      )
    cat(
      paste0(
        'Error rate plot for the ',
        my_direction,
        ' read of primer pair ',
        my_primer_pair_id,
        ' \n'
      )
    )
    plot_errors <-
      dada2::plotErrors(
        error_profile,
        nominalQ = nominalQ,
        obs = obs,
        err_out = err_out,
        err_in = err_in
      )
    ggsave(
      plot_errors,
      filename = 'error_plots.pdf',
      path = directory_path,
      width = 8,
      height = 8
    )
    asv_data <-
      dada2::dada(
        fastq_paths,
        err = error_profile,
        multithread = multithread,
        pool = pool,
        selfConsist = selfConsist
      )

    return(asv_data)

  }

#' Function function to infer ASVs, for multiple loci
#'
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primer_info files
#' 
#' @inheritParams infer_asvs
#' @param denoised_data_path Path to saved intermediate denoised data
#' @keywords internal
infer_asv_command <-
  function(directory_path,
           directory_path_temp,
           multithread = FALSE,
           nbases = 1e+08,
           errorEstimationFunction = loessErrfun,
           randomize = FALSE,
           MAX_CONSIST = 10,
           OMEGA_C = 0,
           qualityType = "Auto",
           nominalQ = FALSE,
           obs = TRUE,
           err_out = TRUE,
           err_in = FALSE,
           pool = FALSE,
           selfConsist = FALSE,
           verbose = FALSE) {
    denoised_data_path <-
      file.path(directory_path_temp, "Denoised_data.Rdata")
    if (file.exists(denoised_data_path)) {
      print("File already exists")
    } else {
      run_dada <- function(direction) {
        lapply(unique(data_tables$cutadapt_data$primer_name), function(primer_name)
          infer_asvs(
            directory_path,
            primer_name,
            direction,
            multithread = multithread,
            nbases = nbases,
            errorEstimationFunction = errorEstimationFunction,
            randomize = randomize,
            MAX_CONSIST = MAX_CONSIST,
            OMEGA_C = OMEGA_C,
            qualityType = qualityType,
            nominalQ = nominalQ,
            obs = obs,
            err_out = err_out,
            err_in = err_in,
            pool = pool,
            selfConsist = selfConsist,
            verbose = verbose
          )) %>%
          unlist(recursive = FALSE)
      }
      dada_forward <- run_dada("Forward")
      dada_reverse <- run_dada("Reverse")
      save(dada_forward, dada_reverse, file = denoised_data_path)
      print("File is now saved")
    }
  }

#' Merge forward and reverse reads
#'
#' @inheritParams dada2::mergePairs
#' @param directory_path A path to the intermediate folder and directory
#' @param merged_read_data_path Path to R data file containing merged read data
#' @return merged_reads Intermediate merged read R data file
#' @keywords internal
merge_reads_command <-
  function(directory_path,
           minOverlap = 12,
           maxMismatch = 0,
           returnRejects = FALSE,
           justConcatenate = FALSE,
           trimOverhang = FALSE,
           verbose = FALSE) {
    denoised_data_path <-
      file.path(directory_path_temp, "Denoised_data.Rdata")
    load(denoised_data_path) #incorporate into function
    merged_read_data_path <-
      file.path(directory_path_temp, "Merged_reads.Rdata")
    formatted_ref_dir <-
      if (file.exists(merged_read_data_path)) {
        load(merged_read_data_path)
        return(merged_reads)
      } else {
        merged_reads <- dada2::mergePairs(
          dadaF = dada_forward,
          derepF = file.path(
            directory_path_temp,
            'filtered_sequences',
            names(dada_forward)
          ),
          dadaR = dada_reverse,
          derepR = file.path(
            directory_path_temp,
            'filtered_sequences',
            names(dada_reverse)
          ),
          minOverlap = minOverlap,
          maxMismatch = maxMismatch,
          returnRejects = returnRejects,
          justConcatenate = justConcatenate,
          verbose = verbose
        )
        names(merged_reads) <- gsub(".fastq.gz", "",
                                    gsub("R1_", "", names(merged_reads), fixed = TRUE))
        save(merged_reads, file = merged_read_data_path)
        return(merged_reads)
      }
  }

#' Count overlap to see how well the reads were merged
#'
#' @param merged_reads Intermediate merged read R data file
#' @return A plot describing how well reads merged and information on overlap between reads
countOverlap <- function(merged_reads, directory_path) {
  non_empty_merged_reads <-
    merged_reads[map_dbl(merged_reads, nrow) > 0]
  non_empty_merged_reads
  merge_data <- non_empty_merged_reads %>%
    bind_rows() %>%
    mutate(fullsample_name = rep(
      names(non_empty_merged_reads),
      map_int(non_empty_merged_reads, nrow)
    )) %>%
    as_tibble() #check best practices

  data_tables$cutadapt_data$fullsample_name <-
    paste0(data_tables$cutadapt_data$file_id,
           '_',
           data_tables$cutadapt_data$primer_name)
  cutadapt_short <- data_tables$cutadapt_data %>%
    select(fullsample_name, sample_name, primer_name)
  merge_data2 <- merge_data %>%
    left_join(cutadapt_short, by = c("fullsample_name" = "fullsample_name"))
  merge_data2 <- mutate(
    merge_data2,
    overlap = nmatch + nmismatch,
    mismatch = nmismatch + nindel,
    identity = (overlap - mismatch) / overlap
  )
  merge_plot <- merge_data2 %>%
    select(primer_name, mismatch, accept, overlap) %>%
    rename(
      'locus' = primer_name,
      'Mismatches and Indels' = mismatch,
      'Merged' = accept,
      'Overlap Length' = overlap
    ) %>%
    gather(key = 'stat', value = 'value',-locus,-Merged) %>%
    ggplot(aes(x = value, fill = Merged)) +
    facet_grid(locus ~ stat, scales = 'free') +
    geom_histogram(bins = 50) +
    scale_fill_viridis_d(begin = 0.8, end = 0.2) +
    labs(x = '', y = 'ASV count', fill = 'Merged') +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  ggsave(
    merge_plot,
    filename = 'read_merging_info.pdf',
    path = directory_path,
    width = 8,
    height = 8
  )
  merge_plot
}

#' Make ASV sequence matrix
#'
#' @inheritParams dada2::makeSequenceTable
#' @param merged_reads Intermediate merged read R data file
#' @return raw_seqtab
#' @keywords internal
makeSeqtab <- function(merged_reads, orderBy = "abundance") {
  raw_seqtab <- makeSequenceTable(merged_reads, orderBy = orderBy)
  return(raw_seqtab)
}

#' Quality filtering to remove chimeras and short sequences
#'
#' @inheritParams dada2::removeBimeraDenovo
#' @param raw_seqtab R data file with raw sequence data prior to removal of chimeras
#' @param min_asv_length The minimum cutoff of ASVs to be retained
#' @return asv_abund_matrix The returned final ASV abundance matrix
#' @keywords internal
make_abund_matrix <-
  function(raw_seqtab,
           method = "consensus",
           verbose = FALSE,
           min_asv_length = 50) {
    seqtab.nochim <-
      dada2::removeBimeraDenovo(raw_seqtab, method = method, verbose = verbose)
    asv_abund_matrix <-
      seqtab.nochim[, nchar(colnames(seqtab.nochim)) >= min_asv_length]
    asvabund_matrix_path <-
      file.path(directory_path_temp, "asvabund_matrixDADA2.Rdata")
    save(asv_abund_matrix, file = asvabund_matrix_path)
    return(asv_abund_matrix)
  }

#' Plots a histogram of read length counts of all sequences within the ASV matrix
#'
#' @param asv_abund_matrix The returned final ASV abundance matrix
#' @keywords internal
#' @return histogram with read length counts of all sequences within ASV matrix
#' @keywords internal
make_seqhist <- function(asv_abund_matrix) {
  ab_matrix <- nchar(getSequences(asv_abund_matrix))
  data1 <- data.frame(ab_matrix)
  hist_plot <-
    ggplot2::qplot(
      ab_matrix,
      data = data1,
      geom = "histogram",
      xlab = 'Length of sequence (bp)',
      ylab = 'Counts',
      main = 'Read length counts of all sequences within the ASV matrix'
    )
  ggsave(
    hist_plot,
    filename = "asv_seqlength_plot.pdf",
    path = directory_path,
    width = 8,
    height = 8
  )
}

