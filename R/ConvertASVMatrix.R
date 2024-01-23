#Add function to combine final matrix if that is desired (need one step to make sure that  no redundant ASVs)

#' Filter ASV abundance matrix and convert to taxmap object
#'
#' @param min_read_depth Threshold for ASVs to remove if number of
#' reads is less than this value across all samples-todo check on this
#' @param minimum_bootstrap Threshold for bootstrap support value
#' for taxonomic assignments. Below designated minimum bootstrap
#' threshold, taxnomoic assignments will be set to N/A.
#' @param pid_species If percent identity is below this value, at the species
#' level, taxonomic assignment will be set to N/A
#' @param pid_genus If percent identity is below this value, at the genus
#' level, taxonomic assignment will be set to N/A
#' @param pid_family If percent identity is below this value, at the family
#' level, taxonomic assignment will be set to N/A
#' @param save_outputs Logical, indicating whether to save the taxmap object. Default is FALSE.
#' @return ASV matrix converted to taxmap object
#' @export
#' @examples 
#' Convert final matrix to taxmap and phyloseq objects for downstream analysis steps
#' prepare_reads(maxN = 0, data_directory = "~/demulticoder/inst/extdata", output_directory = "~/testing_package", tempdir_id = "run1", overwrite_existing = TRUE)
#' cut_trim(analysis_setup,cutadapt_path="/opt/homebrew/bin/cutadapt", overwrite_existing = TRUE)
#' make_asv_abund_matrix(analysis_setup, overwrite_existing = TRUE)
#' assignTax(analysis_setup,asv_abund_matrix, retrieve_files=TRUE, overwrite_existing=TRUE
#' asv_matrix_to_taxmap_phyloseq(save_outputs=TRUE)

asv_matrix_to_taxmap_phyloseq <- function(min_read_depth = 0, minimum_bootstrap = 0, pid_species = 0, pid_genus = 0, pid_family = 0, save_outputs = FALSE) {
  dir_paths <- analysis_setup$dir_paths
  data_tables <- analysis_setup$data_tables
  directory_path <- dir_paths$output_directory
  
  files <- list.files(path = directory_path, pattern = "^final_asv_abundance_matrix_.*\\.csv$", full.names = TRUE)
  suffixes <- gsub("^final_asv_abundance_matrix_(.*)\\.csv$", "\\1", basename(files))
  unique_suffixes <- unique(suffixes)
  
  for (suffix in unique_suffixes) {
    abundance <- read_csv(file.path(directory_path, paste0('final_asv_abundance_matrix_', suffix, '.csv')), show_col_types = FALSE)
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
    
    obj_dada$data$otu_table = obj_dada$data$abund[, -2:-4]
    obj_dada$data$otu_table$otu_id = paste0('ASV', 1:nrow(obj_dada$data$otu_table))
    obj_dada$data$sample_data = data_tables$metadata
    
    phylo_obj <- metacoder::as_phyloseq(
      obj_dada,
      sample_data = obj_dada$data$sample_data,
      sample_id_col = "samplename_barcode"
    )
    
    # Use unique names for saving to avoid overwriting
    taxmap_path <- file.path(directory_path, paste0("taxmap_obj_", suffix, ".RData"))
    phyloseq_path <- file.path(directory_path, paste0("phyloseq_obj_", suffix, ".RData"))
    
    assign(paste0("obj_dada_", suffix), obj_dada, envir = .GlobalEnv)
    assign(paste0("phyloseq_obj_", suffix), phylo_obj, envir = .GlobalEnv)
    
    if (save_outputs == TRUE) {
      save(obj_dada, file = taxmap_path)
      save(phylo_obj, file = phyloseq_path)
      
      cat("For", suffix, "dataset", "\n")
      cat("Taxmap object saved in:", taxmap_path, "\n")
      cat("Phyloseq object saved in:", phyloseq_path, "\n")
      cat("ASVs filtered by minimum read depth:", min_read_depth, "\n")
      cat("For taxonomic assignments, if minimum bootstrap was set to:", minimum_bootstrap, "assignments were set to NA", "\n")
      cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    }
  }
}