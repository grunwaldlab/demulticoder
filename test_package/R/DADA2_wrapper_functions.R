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
#'x
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

#' Retrieve the paths of the filtered and trimmed Fastq ifiles
#'
#' @param my_direction
#' @param my_primer_pair_id
#' @param cutadapt_data Intermediate_data folder with trimmed and filtered reads for each sample
#'
#' @return
#' @export
#'
#' @examples
get_fastq_paths <- function(my_direction, my_primer_pair_id) {
  returnList$cutadapt_data %>%
    filter(direction == my_direction, primer_name == my_primer_pair_id, file.exists(filtered_path)) %>%
    pull(filtered_path)
}

#' Core DADA2 function to learn errors and infer ASVs
#' @inheritParams dada2::learnErrors
#' @inheritParams dada2::dada
#' @inheritParams dada2::plot_errors
#' @param my_primer_pair_id
#' @param my_direction
#'
#' @return
#' @export
#'
#' @examples
infer_asvs <-function(my_primer_pair_id, my_direction, multithread=FALSE,nbases = 1e+08, errorEstimationFunction = loessErrfun, randomize=FALSE, MAX_CONSIST=10, OMEGA_C=0, qualityType="Auto", nominalQ = FALSE, obs=TRUE, err_out=TRUE, err_in=FALSE, pool=FALSE, selfConsist=FALSE, verbose=FALSE){ #may need to adjust the parameters
  fastq_paths <- get_fastq_paths(my_direction, my_primer_pair_id)
  error_profile <- dada2::learnErrors(fastq_paths, multithread = multithread, ,nbases=nbases,errorEstimationFunction= errorEstimationFunction, randomize=randomize, MAX_CONSIST=MAX_CONSIST, OMEGA_C= OMEGA_C, qualityType=qualityType, verbose=verbose)
  #if (error_profile) {
  cat(paste0('Error rate plot for the ', my_direction, ' read of primer pair ', my_primer_pair_id, ' \n'))
  plot_errors<-dada2::plotErrors(error_profile, nominalQ = nominalQ, obs=obs, err_out=err_out, err_in=err_in)
  ggsave(plot_errors, filename = 'error_plots.pdf', path = intermediate_path, width = 8, height = 8)
  asv_data <- dada2::dada(fastq_paths, err = error_profile, multithread = multithread, pool=pool, selfConsist=selfConsist)
  
  return(asv_data)
  
}

#' Function function to infer ASVs, for multiple loci
#' @inheritParams infer_asvs
#' @param intermediate_path A path to the intermediate folder and directory
#' @param cutadapt_data
#' @param denoised_data_path
#'
#' @return
#' @export
#'
#' @examples
infer_asv_command <-function(intermediate_path, cutadapt_data, multithread=FALSE,nbases = 1e+08, errorEstimationFunction = loessErrfun, randomize=FALSE, MAX_CONSIST=10, OMEGA_C=0, qualityType="Auto", nominalQ = FALSE, obs=TRUE, err_out=TRUE, err_in=FALSE, pool=FALSE, selfConsist=FALSE, verbose=FALSE){
  denoised_data_path <- file.path(intermediate_path, "Denoised_data.Rdata")
  if (file.exists(denoised_data_path)) {
    print("File already exists")
  } else {
    run_dada <- function(direction) {
      lapply(unique(returnList$cutadapt_data$primer_name), function(primer_name) infer_asvs(primer_name, direction,multithread = multithread, nbases=nbases,errorEstimationFunction= errorEstimationFunction, randomize=randomize, MAX_CONSIST=MAX_CONSIST, OMEGA_C= OMEGA_C, qualityType=qualityType, nominalQ = nominalQ, obs=obs, err_out=err_out, err_in=err_in,pool=pool, selfConsist=selfConsist, verbose=verbose)) %>%
        unlist(recursive = FALSE)
    }
    dada_forward <- run_dada("Forward")
    dada_reverse <- run_dada("Reverse")
    save(dada_forward, dada_reverse, file = denoised_data_path)
    print("File is now saved")
  }
}

#' Merge forward and reverse reads ###CHECK
#' @inheritParams dada2::mergePairs
#' @param intermediate_path A path to the intermediate folder and directory
#' @param merged_read_data_path
#'
#' @return
#' @export
#'
#' @examples
merge_reads_command <- function(intermediate_path, minOverlap=12, maxMismatch=0, returnRejects=FALSE, justConcatenate=FALSE, trimOverhang=FALSE, verbose=FALSE){
  denoised_data_path <- file.path(intermediate_path, "Denoised_data.Rdata")
  load(denoised_data_path) #incorporate into function
  merged_read_data_path <- file.path(intermediate_path, "Merged_reads.Rdata")
  formatted_ref_dir <-
    if (file.exists(merged_read_data_path)) {
      load(merged_read_data_path)
      return(merged_reads)
    } else {
      merged_reads <- dada2::mergePairs(dadaF = dada_forward,
                                        derepF = file.path(intermediate_path, 'filtered_sequences', names(dada_forward)),
                                        dadaR = dada_reverse,
                                        derepR = file.path(intermediate_path, 'filtered_sequences', names(dada_reverse)),
                                        minOverlap = minOverlap,
                                        maxMismatch = maxMismatch,
                                        returnRejects = returnRejects,
                                        justConcatenate = justConcatenate,
                                        verbose = verbose)
      names(merged_reads) <- gsub(".fastq.gz", "", names(merged_reads), fixed = TRUE)
      save(merged_reads, file = merged_read_data_path)
      return(merged_reads)
    }
}

#' Count overlap to see how well the reads were merged
#'
#' @param merged_reads
#'
#' @return
#' @export
#'
#' @examples
countOverlap <- function(merged_reads){
  non_empty_merged_reads <- merged_reads[map_dbl(merged_reads, nrow) > 0]
  non_empty_merged_reads
  merge_data <- non_empty_merged_reads %>%
    bind_rows() %>%
    mutate(fullsample_id = rep(names(non_empty_merged_reads), map_int(non_empty_merged_reads, nrow))) %>%
    as_tibble() #check best practices
  
  returnList$cutadapt_data$fullsample_id <- paste0(returnList$cutadapt_data$file_id,'_', returnList$cutadapt_data$primer_name)
  cutadapt_short <- returnList$cutadapt_data %>%
    select(fullsample_id, sample_id, primer_name)
  merge_data2 <- merge_data %>%
    left_join(cutadapt_short, by = c("fullsample_id" = "fullsample_id"))
  merge_data2 <- merge_data %>%
    left_join(cutadapt_short, by = c("fullsample_id" = "fullsample_id"))
  merge_data2 <- mutate(merge_data2,
                        overlap = nmatch + nmismatch,
                        mismatch = nmismatch + nindel,
                        identity = (overlap - mismatch) / overlap)
  merge_plot <- merge_data2 %>%
    select(primer_name, mismatch, accept, overlap) %>%
    rename('locus' = primer_name, 'Mismatches and Indels' = mismatch, 'Merged' = accept, 'Overlap Length' = overlap) %>%
    gather(key = 'stat', value = 'value', -locus, -Merged) %>%
    ggplot(aes(x = value, fill = Merged)) +
    facet_grid(locus ~ stat, scales = 'free') +
    geom_histogram(bins = 50) +
    scale_fill_viridis_d(begin = 0.8, end = 0.2) +
    labs(x = '', y = 'ASV count', fill = 'Merged') +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position="bottom")
  ggsave(merge_plot, filename = 'read_merging.jpg', path = intermediate_path, width = 8, height = 8)
  merge_plot
}


#' Make ASV sequence table
#' @inheritParams dada2::makeSequenceTable
#' @param merged_reads
#'
#' @return
#' @export
#'
#' @examples
makeSeqtab<-function(merged_reads, orderBy="abundance"){
  raw_seqtab<-makeSequenceTable(merged_reads, orderBy=orderBy)
  return(raw_seqtab)
}

#Remove chimeras and short sequences
#' Quality filtering to remove chimeras and short sequences
#' @inheritParams dada2::removeBimeraDenovo
#' @param seqtab
#' @param min_asv_length
#'
#' @return
#' @export
#'
#' @examples
make_abund_table<-function(raw_seqtab, method="consensus", verbose=FALSE, min_asv_length=50){
  seqtab.nochim <- dada2::removeBimeraDenovo(raw_seqtab, method=method, verbose=verbose)
  asv_abund_table <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) >= min_asv_length]
  asvabund_table_path <- file.path(intermediate_path, "asvabund_tableDADA2.Rdata")
  save(asv_abund_table, file = asvabund_table_path)
  return(asv_abund_table)
}

#why  not showing up? 
#' Plots a histogram of read length counts of all sequences within the ASV table
#'
#' @param asv_abund_table
#'
#' @return
#' @export
#'
#' @examples
make_seqhist <- function(asv_abund_table){
  ab_table<-nchar(getSequences(asv_abund_table))
  data1<-data.frame(ab_table)
  hist_plot<-ggplot2::qplot(ab_table, data=data1, geom="histogram", xlab = 'Length of sequence (bp)', ylab='Counts', main='Read length counts of all sequences within the ASV table')
  ggsave(hist_plot, filename = "asv_hist_plot.pdf", path = intermediate_path, width = 8, height = 8)
}

#' Prepare final ASV abundance matrix
#'
#' @param Intermediate_data folder with trimmed and filtered reads for each sample
#' @param asv_abund_table
#' @param locus
#'
#' @return
#' @export
#'
#' @examples
prep_abund_table <-function(cutadapt_data, asv_abund_table, locus){
  #rownames(asv_abund_table) <- sub(rownames(asv_abund_table), pattern = ".fastq.gz", replacement = "")
  returnList$cutadapt_data$file_id_primer <- paste0(returnList$cutadapt_data$file_id, "_", returnList$cutadapt_data$primer_name)
  asv_abund_matrix<- asv_abund_table[rownames(asv_abund_table) %in% returnList$cutadapt_data$file_id_primer[returnList$cutadapt_data$primer_name == locus], ]
  return(asv_abund_matrix)
}


#' Still needs improvement. Could solve problem by merging fungal and oomycete databases for assign taxonomy steps. In order to assign taxonomy for separate loci, need to separate abundance tables first
#'
#' @param abund_asv_its
#' @param abund_asv_rps10
#' @param intermediate_path A path to the intermediate folder and directory
#' @param asv_abund_matrix
#'
#' @return
#' @export
#'
#' @examples
separate_abund_table <- function(abund_asv_its, abund_asv_rps10, intermediate_path, asv_abund_matrix){
  separate_abund_path <- file.path(intermediate_path, "Separate_abund.Rdata")
  in_both <- colSums(abund_asv_its) != 0 & colSums(abund_asv_rps10) != 0
  assign_to_its <- in_both & colSums(abund_asv_rps10) > colSums(abund_asv_its)
  assign_to_rps10 <- in_both & colSums(abund_asv_rps10) < colSums(abund_asv_its)
  is_its <- (colSums(abund_asv_its) != 0 & colSums(abund_asv_rps10) == 0) | assign_to_its
  is_rps10 <- (colSums(abund_asv_rps10) != 0 & colSums(abund_asv_its) == 0) | assign_to_rps10
  abund_asv_its <- abund_asv_its[ , is_its]
  abund_asv_rps10 <- abund_asv_rps10[ , is_rps10]
  #The number of ASVs left in the two groups should sum to the total number of ASVs, since there should be no overlap.
  save(abund_asv_its, abund_asv_rps10, file = separate_abund_path)
  stopifnot(ncol(abund_asv_its) + ncol(abund_asv_rps10) == ncol(asv_abund_matrix)) #make more messaging
}

#This function is unwieldy if not enough 
#' Assign taxonomy
#' @inheritParams dada2::assignTaxonomy
#' @param abund_asv_table
#' @param ref_database
#' @param taxresults_file
#'
#' @return
#' @export
#'
#' @examples
assign_taxonomyDada2<-function(abund_asv_table, ref_database, taxresults_file, minBoot=0, tryRC=FALSE, verbose=FALSE, multithread=TRUE){
  tax_results<- dada2::assignTaxonomy(abund_asv_table,
                                      refFasta = file.path(intermediate_path, 'reference_databases', ref_database),
                                      taxLevels = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Reference"),
                                      minBoot = minBoot,
                                      tryRC = tryRC,
                                      outputBootstraps = TRUE,
                                      multithread = multithread)
  tax_table_path <- file.path(intermediate_path, taxresults_file)
  save(tax_results, file = tax_table_path)
  return(tax_results)
  #save these results
}
#Put these functions in separate directory. 
#' Align ASV sequences to reference sequences from database to get percent ID. STart by retrieving reference sequences.
#'
#' @param tax_results
#' @param db
#'
#' @return
#' @export
#'
#' @examples
get_ref_seq <- function(tax_results, db) {
  ref_i <- as.integer(str_match(tax_results$tax[, 'Reference'], '^.+_([0-9]+)$')[ ,2])
  db[ref_i]
}

#' Align ASV sequences to reference sequences from database to get percent ID. Align sequences to reference sequences from database.
#'
#' @param ref
#' @param asv
#'
#' @return
#' @export
#'
#' @examples
get_align_pid <- function(ref, asv) {
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  align <-  pairwiseAlignment(pattern = asv, subject = ref, type = 'global-local')
  is_match <- strsplit(as.character(align@pattern), '')[[1]] == strsplit(as.character(align@subject), '')[[1]]
  sum(is_match) / length(is_match)
}

#' Align ASV sequences to reference sequences from database to get percent ID. Get percent identities.
#'
#' @param tax_results
#'
#' @return
#' @export
#'
#' @examples
get_pids <- function(tax_results, db) {
  db_seqs <- read_fasta(file.path(intermediate_path, "reference_databases", db))
  ref_seqs <- get_ref_seq(tax_results, db_seqs)
  asv_seqs <- rownames(tax_results$tax)
  asv_ref_align_normal <- map2(asv_seqs, ref_seqs, pairwiseAlignment, type = 'global-local')
  asv_pids_normal <- map_dbl(asv_ref_align_normal, pid, type = "PID1")
  asv_ref_align_revcomp <- map2(asv_seqs, rev_comp(ref_seqs), pairwiseAlignment, type = 'global-local')
  asv_pids_revcomp <- map_dbl(asv_ref_align_revcomp, pid, type = "PID1")
  asv_ref_align <- ifelse(asv_pids_normal >= asv_pids_revcomp, asv_ref_align_normal, asv_ref_align_revcomp)
  asv_pids <- ifelse(asv_pids_normal >= asv_pids_revcomp, asv_pids_normal, asv_pids_revcomp)
  paste0(map_chr(seq_along(asv_ref_align), function(i) {
    paste0('asv: ',  asv_seqs[[i]], '\n',
           'ref: ', names(ref_seqs)[[i]], '\n',
           'pid: ', asv_pids[[i]], '\n',
           asv_ref_align[[i]]@pattern, '\n',
           asv_ref_align[[i]]@subject, '\n')
  }), collapse = '==================================================\n') %>%
    write_lines(file = file.path(intermediate_path, 'dada2_asv_alignments.txt'))
  return(asv_pids)
}

#' Add PID and bootstrap values to tax result.
#'
#' @param tax_results
#' @param pid
#'
#' @return
#' @export
#'
#' @examples
add_pid_to_tax <- function(tax_results, pid) {
  tax_results$tax <- cbind(tax_results$tax, ASV = rownames(tax_results$tax))
  tax_results$boot <- cbind(tax_results$boot, ASV = pid)
  return(tax_results)
}

#' Combine taxonomic assignments and bootstrap values for each locus into single classfication vector
#'
#' @param tax_results
#'
#' @return
#' @export
#'
#' @examples
assignTax_as_char <- function(tax_results) {
  tax_table_path <- file.path(intermediate_path, "Final_tax_table.Rdata")
  tax_out <- vapply(1:nrow(tax_results$tax), FUN.VALUE = character(1), function(i) {
    paste(tax_results$tax[i, ],
          tax_results$boot[i, ],
          colnames(tax_results$tax),
          sep = '--', collapse = ';')
  })
  names(tax_out) <- rownames(tax_results$tax)
  save(tax_out, file = tax_table_path)
  return(tax_out)
  #check
  stopifnot(all(names(tax_out) %in% colnames(asv_abund_matrix)))
  stopifnot(all(! duplicated(names(tax_out))))
}

#' Format ASV abundance table
#'
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#' @param seq_tax_asv An amplified sequence variants table with taxonomic information
#'
#' @return
#' @export
#'
#' @examples
#Reformat ASV table
format_abund_matrix <- function(asv_abund_matrix, seq_tax_asv) {
  formatted_abund_asv <- t(asv_abund_matrix)
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



#' Final inventory of read counts after each step from input to removal of chimeras. This function deals with if you have more than one sample. TODO optimize for one sample
#'
#' @param asv_abund_matrix An abundance matrix containing amplified sequence variants
#'
#' 
#' @return
#' @export
#'
#' @examples
dada2_readcounts_multi_sample <- function(asv_abund_matrix) {
  filter_results_path<-file.path(intermediate_path, "filter_results.RData")
  load(filter_results_path) #incorporate into function
  denoised_data_path <- file.path(intermediate_path, "Denoised_data.Rdata")
  load(denoised_data_path) #incorporate into function
  merged_read_data_path <- file.path(intermediate_path, "Merged_reads.Rdata")
  load(merged_read_data_path)
  getN <- function(x) sum(dada2::getUniques(x))
  track <- cbind(filter_results, sapply(dada_forward, getN), sapply(dada_reverse, getN), sapply(merged_reads, getN), rowSums(asv_abund_matrix))
  track <- cbind(rownames(track), data.frame(track, row.names=NULL))
  colnames(track) <- c("sample_name", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  track_read_counts_path <- file.path(intermediate_path, "track_reads.csv")
  write_csv(track, track_read_counts_path)
  print(track)
}