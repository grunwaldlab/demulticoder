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

#' Retrieve the paths of the filtered and trimmed Fastq files
#'
#' @param my_direction Whether primer is in forward or reverse direction
#' @param my_primer_pair_id The the specific barcode id 
#' @param cutadapt_data directory_data folder with trimmed and filtered reads for each sample
#' @keywords internal
get_fastq_paths <- function(my_direction, my_primer_pair_id) {
  filtered_paths <- character()
  
  for (i in seq_along(data_tables$cutadapt_data$direction)) {
    if (data_tables$cutadapt_data$direction[i] == my_direction &&
        data_tables$cutadapt_data$primer_name[i] == my_primer_pair_id &&
        file.exists(data_tables$cutadapt_data$filtered_path[i])) {
      filtered_paths <- c(filtered_paths, data_tables$cutadapt_data$filtered_path[i])
    }
  }
  
  filtered_paths
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
        dada_output <- lapply(unique(data_tables$cutadapt_data$primer_name), function(primer_name)
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
          ))
        unlist(dada_output, recursive = FALSE)
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
#' Count overlap to see how well the reads were merged
#'
#' @param merged_reads Intermediate merged read R data file
#' @param directory_path Directory path to save the plot
#' @return A plot describing how well reads merged and information on overlap between reads
countOverlap <- function(merged_reads, directory_path) {
  library(ggplot2)
  
  non_empty_merged_reads <- merged_reads[sapply(merged_reads, nrow) > 0]
  
  merge_data <- do.call(rbind, non_empty_merged_reads)
  merge_data$samplename_barcode <- rep(names(non_empty_merged_reads), sapply(non_empty_merged_reads, nrow))
  
  merge_data2 <- merge_data[merge_data$samplename_barcode %in% data_tables$cutadapt_data$samplename_barcode, ]
  merge_data2 <- merge(merge_data2, data_tables$cutadapt_data, by = "samplename_barcode")
  
  merge_data2$overlap <- merge_data2$nmatch + merge_data2$nmismatch
  merge_data2$mismatch <- merge_data2$nmismatch + merge_data2$nindel
  merge_data2$identity <- (merge_data2$overlap - merge_data2$mismatch) / merge_data2$overlap
  
  # Make dataframe long
  merge_plot <- data.frame(
    locus = merge_data2$primer_name,
    Mismatches_and_Indels = merge_data2$mismatch,
    Merged = merge_data2$accept,
    Overlap_Length = merge_data2$overlap
  )
  
  long_df <- stack(merge_plot[,c("Mismatches_and_Indels", "Overlap_Length")])
  merge_plot <- cbind(merge_plot[,"locus"], merge_plot[,"Merged"], long_df)
  names(merge_plot) <- c("locus", "Merged", "value", "stat")
  
  merge_plot_output <- ggplot(merge_plot, aes(x = value, fill = Merged)) +
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
    merge_plot_output,
    filename = 'read_merging_info.pdf',
    path = directory_path,
    width = 8,
    height = 8
  )
  
  merge_plot_output
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



#Make separate plots for each barcode? 
#' Plots a histogram of read length counts of all sequences within the ASV matrix
#'
#' @param asv_abund_matrix The returned final ASV abundance matrix
#' @keywords internal
#' @return histogram with read length counts of all sequences within ASV matrix
#' @keywords internal
make_seqhist <- function(asv_abund_matrix) {
  matrix_row <- unique(gsub(".*_", "", rownames(asv_abund_matrix)))
  
  for (locus in matrix_row) {
    locus_ab_matrix <- nchar(getSequences(asv_abund_matrix[, grepl(paste0("_", locus), rownames(asv_abund_matrix))]))
    data <- data.frame(locus_ab_matrix)
    
    hist_plot <- ggplot2::qplot(
      locus_ab_matrix,
      data = data,
      geom = "histogram",
      xlab = 'Length of sequence (bp)',
      ylab = 'Counts',
      main = paste("Read length counts of ASVs produced from", locus, "amplicon data")
    )
    
    ggsave(
      hist_plot,
      filename = paste("asv_seqlength_plot_",locus,".pdf", sep = ""),
      path = directory_path,
      width = 8,
      height = 8
    )
  }
}
