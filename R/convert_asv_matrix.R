#' Filter ASV abundance matrix and convert to taxmap object
#'
#' @param analysis_setup analysis_setup An object containing directory paths and
#'   data tables, produced by the `prepare_reads` function
#' @param min_read_depth ASV filter parameter. If mean read depth of across all
#'   samples is less than this threshold, ASV will be filtered. reads is less
#'   than this value across all samples
#' @param minimum_bootstrap Threshold for bootstrap support value for taxonomic
#'   assignments. Below designated minimum bootstrap threshold, taxnomoic
#'   assignments will be set to N/A
#' @param save_outputs Logical, indicating whether to save the taxmap object.
#'   Default is FALSE.
#' @return ASV matrix converted to taxmap object
#' @export
#' @examples
#' Convert final matrix to taxmap and phyloseq objects for downstream analysis steps
#' prepare_reads(data_directory = "inst/extdata", output_directory = "test_data", tempdir_id = "demulticoder_run", overwrite_existing = FALSE)
#' cut_trim(analysis_setup,cutadapt_path="/opt/homebrew/bin/cutadapt", overwrite_existing = FALSE)
#' make_asv_abund_matrix(analysis_setup, overwrite_existing = FALSE)
#' assign_tax(analysis_setup,asv_abund_matrix, retrieve_files=FALSE, overwrite_existing=FALSE)
#' convert_asv_matrix_to_objs(analysis_setup, save_outputs=FALSE)
convert_asv_matrix_to_objs <- function(analysis_setup, min_read_depth = 0, minimum_bootstrap = 0, save_outputs = FALSE, overwrite = FALSE) {
  data_tables <- analysis_setup$data_tables
  output_directory_path <- analysis_setup$directory_paths$output_directory

  files <- list.files(path = output_directory_path, pattern = "^final_asv_abundance_matrix_.*\\.csv$", full.names = TRUE)
  suffixes <- gsub("^final_asv_abundance_matrix_(.*)\\.csv$", "\\1", basename(files))
  unique_suffixes <- unique(suffixes)
  
  for (suffix in unique_suffixes) {
    taxmap_name <- paste0("obj_dada_", suffix)
    phyloseq_name <- paste0("phylo_obj_", suffix)
    
    taxmap_path <- file.path(output_directory_path, paste0(taxmap_name, ".RData"))
    phyloseq_path <- file.path(output_directory_path, paste0(phyloseq_name, ".RData"))
    
    if (file.exists(taxmap_path) && file.exists(phyloseq_path) && !overwrite) {
      cat("For", suffix, "dataset", "\n")
      cat("Files already exist:", taxmap_path, "and", phyloseq_path, "\n")
      cat("To overwrite, set overwrite = TRUE\n")
      cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    } else {
      abundance <- read_csv(file.path(output_directory_path, paste0('final_asv_abundance_matrix_', suffix, '.csv')))
      is_low_abund <- rowSums(abundance[, grepl(paste0("_", suffix, "$"), colnames(abundance))]) < min_read_depth
      abundance <- filter(abundance, !is_low_abund)
      abundance$dada2_tax <- map_chr(strsplit(abundance$dada2_tax, ';'), function(x) {
        map_chr(strsplit(x, '--'), function(parts) {
          if (parts[[3]] != "ASV" && as.numeric(parts[[2]]) <= minimum_bootstrap) {
            parts[[1]] <- "Unknown"
          }
          return(paste(parts, collapse = '--'))
        }) %>% paste(collapse = ';')
      })
      obj_dada <- metacoder::parse_tax_data(abundance, class_cols = 'dada2_tax', class_sep = ';', include_tax_data = TRUE,
                                            class_regex = '^(.+)--(.+)--(.+)$',
                                            class_key = c(taxon = 'taxon_name', boot = 'info', rank = 'taxon_rank'))
      names(obj_dada$data) <- c('abund', 'score')
      
      obj_dada$data$otu_table = obj_dada$data$abund[, -3:-4]
      #obj_dada$data$otu_table$otu_id = paste0('ASV', 1:nrow(obj_dada$data$otu_table))
      obj_dada$data$sample_data = data_tables$metadata
      
      phylo_obj <- metacoder::as_phyloseq(
        obj_dada,
        sample_data = obj_dada$data$sample_data,
        sample_id_col = "samplename_barcode", 
        otu_id_col = "asv_id"
      )
      
      assign(taxmap_name, obj_dada, envir = .GlobalEnv)
      assign(phyloseq_name, phylo_obj, envir = .GlobalEnv)
      
      if (save_outputs == TRUE) {
        save(obj_dada, file = taxmap_path)
        save(phylo_obj, file = phyloseq_path)
        
        cat("For", suffix, "dataset", "\n")
        cat("Taxmap object saved in:", taxmap_path, "\n")
        cat("Phyloseq object saved in:", phyloseq_path, "\n")
        cat("ASVs filtered by minimum read depth:", min_read_depth, "\n")
        cat("For taxonomic assignments, if minimum bootstrap was set to:", minimum_bootstrap, "assignments were set to NA", "\n")
        cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
      } else {
        cat("For", suffix, "dataset", "\n")
        cat("ASV matrix found, but save_outputs is FALSE. Rerun previous part of the pipeline.\n")
        cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
      }
    }
  }
}